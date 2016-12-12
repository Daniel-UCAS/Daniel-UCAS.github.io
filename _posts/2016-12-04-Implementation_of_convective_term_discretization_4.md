---
layout: post
title: Implementation of convective term discretization (4)
author: Daniel-UCAS
date: 2016-12-04
category: OpenFOAM
tags:
- OpenFOAM
- Convetive term
- Discretiazation
excerpt: TVD-NVD framework in OpenFOAM
---

* TOC
{:toc}

### Implementation of limited schemes in OpenFOAM
`OpenFOAM`中限制器格式的实现与常规`CFD`教材(如`HK Versteeg, Malalsekera`的`An Introduction to Computational Fluid Dynamics`第二版, P166-167)中的写法有一点区别，说明如下。
通常`CFD`教材中介绍`TVD`格式的思路是将二阶中心格式写成一阶迎风格式加上一个修正项，然后在修正项前面引入限制函数，即
$$
(\phi_f)_{CD} = \frac{1}{2}(\phi_C + \phi_D) = \phi_C + \frac{1}{2}(\phi_D - \phi_C)  \qquad  (1)
$$
引入限制函数$\psi(r)$后格式写成
$$
(\phi_f)_{TVD} = \phi_C + \frac{1}{2}\psi(r)(\phi_D - \phi_C)  \qquad  (2)
$$
上述中心格式针对的是均匀网格，与`OpenFOAM`中不太一样。`OpenFOAM`中对流项格式离散的基类在[`surfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.C)中(4.x之后版本)，其中最主要的`interpolate`函数的实现为
```cpp
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    return dotInterpolate(geometricOneField(), vf, tlambdas);
}
```
其中调用了`dotInterpolate`函数，相应的主要代码为
```cpp
for (label fi=0; fi<P.size(); fi++)
{
    sfi[fi] = Sfi[fi] & (lambda[fi]*(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]]);
}
```
忽略与`Sfi[fi]`的内积，对应的插值公式为
$$
(\phi_f)_{CD} = \lambda*\phi_P + (1-\lambda)*\phi_N  \qquad  (3)
$$
其中，$\lambda$是插值权重。老版本`OpenFOAM`的`interpolate`函数的实现在[`surfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-3.0.x/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.C)中

```cpp
for (label fi=0; fi<P.size(); fi++)
for (label fi=0; fi<P.size(); fi++)
{
    sfi[fi] = lambda[fi] *(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]];
}
```
直接就是公式(3)。
在显式修正时，[`surfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-3.0.x/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.C)的第三个`interpolate`函数即调用了上述两个入口参数的`interpolate`函数
```cpp
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    ...
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
        = interpolate(vf, weights(vf));
    if (corrected())
    {
        tsf() += correction(vf);
    }
    return tsf;
}
```
其中的`weights(vf)`函数是[`surfaceInterpolationScheme.H`](https://github.com/OpenFOAM/OpenFOAM-3.0.x/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H)中的纯虚函数，对于非限制器格式，每种格式均有一个`weights`函数的实现，对于限制器格式，该虚函数的实现位于[`limitedSurfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.C)中
```cpp
template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::limitedSurfaceInterpolationScheme<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    return this->weights
    (
        phi,
        this->mesh().surfaceInterpolation::weights(),
        this->limiter(phi)
    );
}
```
这一实现首先调用三个入口参数的`weights(phi,this->mesh().surfaceInterpolation::weights(),this->limiter(phi))`函数，该函数又调用了`surfaceInterpolation::weights()`，用以计算中心差分格式的权重，还调用了用于计算限制器$\psi(r)$的限制函数`limiter(phi)`。**需要注意的是，每个不同的限制器格式都会有一个`limiter(phi)`函数，这也是不同限制器格式的主要区别**
先看`surfaceInterpolation::weights()`，程序实现位于[`surfaceInterpolation.C`](https://github.com/OpenFOAM/OpenFOAM-3.0.x/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolation.C)中
```cpp
const Foam::surfaceScalarField&
Foam::surfaceInterpolation::weights() const
{
    if (!weights_)
    {
        makeWeights();
    }

    return (*weights_);
}
...
void Foam::surfaceInterpolation::makeWeights() const
{
    ...
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Sf = mesh_.faceAreas();

    // ... and reference to the internal field of the weighting factors
    scalarField& w = weights.internalField();

    forAll(owner, facei)
    {
        scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
        scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));
        w[facei] = SfdNei/(SfdOwn + SfdNei);
    }

    forAll(mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].makeWeights
        (
            weights.boundaryField()[patchi]
        );
    }
    ...
}
```
注意最后一个`forAll`循环是边界条件的实现，其中调用的`makeWeights`函数在[`fvPatch.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/fba2c4f29b3fa9d58351f18e19ac14679f0703bd/src/finiteVolume/fvMesh/fvPatches/fvPatch/fvPatch.C)中，返回的权重$w = 1.0$。上述被`weights`函数调用的`makeWeights()`函数(倒数第二个`forAll`循环)对应的公式为
$$
\omega_{CD} = \frac{\vec{S}_f \cdot \vec{fN}}{\vec{S}_f \cdot \vec{PN}}
$$
得到的结果是二阶中心差分格式的插值权重。上述公式和下面的示意图均来自`Hrvoje Jasak`的博士论文`P81-82`
![`face interpolation`](https://github.com/Daniel-UCAS/images_for_markdown/blob/master/2016-12-05_face_interpolation.jpg?raw=true)
再看三个入口参数的`weights`函数，其实现也位于[`limitedSurfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.C)中
```cpp
template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::limitedSurfaceInterpolationScheme<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const surfaceScalarField& CDweights,
    tmp<surfaceScalarField> tLimiter
) const
{
    // Note that here the weights field is initialised as the limiter
    // from which the weight is calculated using the limiter value
    surfaceScalarField& Weights = tLimiter.ref();

    scalarField& pWeights = Weights.primitiveFieldRef();

    forAll(pWeights, face)
    {
        pWeights[face] =
            pWeights[face]*CDweights[face]
          + (1.0 - pWeights[face])*pos(faceFlux_[face]);
    }

    surfaceScalarField::Boundary& bWeights =
        Weights.boundaryFieldRef();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];

        const scalarField& pCDweights = CDweights.boundaryField()[patchi];
        const scalarField& pFaceFlux = faceFlux_.boundaryField()[patchi];

        forAll(pWeights, face)
        {
            pWeights[face] =
                pWeights[face]*pCDweights[face]
              + (1.0 - pWeights[face])*pos(pFaceFlux[face]);
        }
    }

    return tLimiter;
}
```
这一函数对应的计算公式为
$$
\lambda = \psi(r) * \omega_{CD} + (1-\psi(r)) \, pos(flux) \qquad  (4)
$$
注意$pos(flux)$函数用于判断流动的方向
$$
pos(x) =\begin{cases}1 & x \geq 0
        \\0 & x \lt 0\end{cases}
$$
考虑$flux \geq 0$的情况，
$$
\lambda = \psi(r) * \omega_{CD} + (1-\psi(r))
        = 1 - \psi(r)(1- \omega_{CD})  \qquad  (5)
$$
将(5)式中的$\lambda$代入(3)式中可以得到`TVD`形式的插值
$$
\begin{aligned}
(\phi_f)_{limited} &= (1 - \psi(r)(1- \omega_{CD}))\phi_P+ \psi(r)(1- \omega_{CD})\phi_N\\
 &= \phi_P + \psi(r)(1-\omega_{CD})(\phi_N - \phi_P) \qquad  (6)
\end{aligned}
$$
综合(2)式和(6)式可知常规`CFD`教材中限制器的写法为
$$
(\phi_f){TVD} = \phi_C + \frac{1}{2}\psi(r)(\phi_D - \phi_C)
$$
`OpenFOAM`中限制器函数的实现为
$$
(\phi_f)_{TVD} = \phi_P + \psi(r)(1-\omega_{CD})(\phi_N - \phi_P)
$$
对于均匀网格，$\omega_{CD}=\frac{1}{2}$，`OpenFOAM`中的程序实现退化为常规`CFD`教材的形式。
备注：`OpenFOAM`中`TVD`格式的实现方式在`Hrvoje Jasak`的博士论文[Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows](http://powerlab.fsb.hr/ped/kturbo/OpenFOAM/docs/HrvojeJasakPhD.pdf)中有较好的论述，现推导如下。
首先`TVD`格式的插值形式为公式`(3.66)`，见于`P98`
$$
(\phi_f)_{TVD} = (\phi_f)_{UD} + \psi(r) ((\phi_f)_{HO} - (\phi_f)_{UD})
$$
其中$\phi_{HO}$是以二阶中心差分格式实现的，见`P82`的`公式(3.19)`
$$
(\phi_f)_{HO} = (\phi_f)_{CD} = f_x \phi_P + (1 - f_x) \phi_N = \omega_{CD} \; \phi_P + (1 - \omega_{CD}) \phi_N
$$
$(\phi_f)_{UD}$是一阶迎风格式，见`P83`的`公式(3.21)`
$$
(\phi_f)_{UD} =\begin{cases}\phi_P & flux \geq 0
        \\ \phi_N & flux \lt 0\end{cases}
$$
同样考虑$flux \geq 0$的情况，将$(\phi_f)_{UD}$和$(\phi_f)_{HO}$代入到$(\phi_f)_{TVD}$中可以得到
$$
\begin{aligned}
(\phi_f)_{TVD} &= \phi_P + \psi(r) (\omega_{CD} \; \phi_P + (1 - \omega_{CD}) \phi_N - \phi_P)\\
 &= \phi_P + \psi(r)(1-\omega_{CD})(\phi_N - \phi_P)
\end{aligned}
$$

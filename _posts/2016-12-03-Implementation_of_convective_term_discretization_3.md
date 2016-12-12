---
layout: post
title: Implementation of convective term discretization (3)
author: Daniel-UCAS
date: 2016-12-03
category: OpenFOAM
tags:
- OpenFOAM
- Convetive term
- Discretiazation
excerpt: Interpretation for interpolation in OpenFOAM.
---

* TOC
{:toc}

### `interpolate` function
`OpenFOAM`中对流项的离散核心问题是利用当前单元P(C)和相邻单元N(下游单元D)的变量值来计算两个单元工友界面(e)上的变量值，这一计算通过插值(interpolate)来实现，称之为重构(reconstruct)。**需要注意的是`OpenFOAM`以非结构方式存储网格，因此一般只能利用两个相邻单元的变量值。**
用于插值的`interpolate`函数位于[`surfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.C)文件中，共有四个版本。前两个是给定单元的权重计算面上的变量值，第一个是给定两个单元(P、N)的权重，第二个是给定一个单元(当前单元P)的权重；另外还有两个用于变量的显式修正，来实现高阶格式。
使用P、N两个节点权重加权的`interpolate`函数计算原理为
$$
\phi_e = \lambda \phi_P + \nu \phi_N
$$
相应的程序实现为
```cpp
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas,
    const tmp<surfaceScalarField>& tys
)
{
    ...
    for (label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = lambda[fi]*vfi[P[fi]] + y[fi]*vfi[N[fi]];
    }

    // Interpolate across coupled patches using given lambdas and ys
    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const fvsPatchScalarField& pY = ys.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            sf.boundaryField()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
              + pY*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }
    ...
}
```

前一个`for`循环计算内部单元面，后一个`forAll`循环处理边界面，`coupled`模块处理并行。
使用当前单元的权函数插值的`interpolate`函数计算原理为
$$
\phi_e = \lambda \phi_P + (1-\lambda) \phi_N = \lambda (\phi_P-\phi_N) + \phi_N
$$
对应的程序实现为
```cpp
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    ...

    for (label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = lambda[fi]*(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]];
    }
    // Interpolate across coupled patches using given lambdas
    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            tsf().boundaryField()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
             + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }
    ...
}
```
增加显式修正的两个`interpolate`函数计算原理相同，都是在插值得到的变量上增加某种修正，唯一的区别是函数入口参数类型，一个是`const GeometricField<Type, fvPatchField, volMesh>&`，另一个是`const tmp<GeometricField<Type, fvPatchField, volMesh> >& `。前一个函数的实现程序为
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
该函数调用了上面的第二个版本的`interpolate`函数进行插值。上面给出的代码为旧版本的程序，2016年以后的`OpenFOAM-4.x`版本中新引入了两个`dotInterpolate`函数，可以避免在计算时额外增加临时变量，以减少在计算速度梯度的散度时所需的内存开销，有时也可以减少计算时间，参考[Added dotInterpolate member-function](http://caefn.com/openfoam/solvers-recent-changes)。

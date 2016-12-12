---
layout: post
title: Implementation of convective term discretization (2)
author: Daniel-UCAS
date: 2016-12-02
category: OpenFOAM
tags:
- OpenFOAM
- Convetive term
- Discretiazation
excerpt: Code implementation of convective term in OpenFOAM.
---

* TOC
{:toc}

### Implementation of non-limited schemes

#### 离散方法

一阶迎风格式示意图如下

![upwind profile](https://github.com/Daniel-UCAS/images_for_markdown/blob/master/2016-11-28_upwind_scheme_profile.jpg?raw=true)


$$
\phi_e =\begin{cases}
        \phi_C & if\, \dot{m}_e>0\\
        \phi_E & if\, \dot{m}_e<0
        \end{cases} \qquad
\phi_w =\begin{cases}
        \phi_C & if\, \dot{m}_e>0\\
        \phi_W & if\, \dot{m}_e<0
        \end{cases}        
$$


#### 程序实现

继承自`limitedSurfaceInterpolationScheme<Type>`类，主要重写了`weights()`函数

```cpp
tmp<surfaceScalarField> weights() const
{
    return pos(this->faceFlux_);
}

```

### 2nd order central difference scheme

#### 离散方

二阶中心差分格式示意图如下

![CD profile](https://github.com/Daniel-UCAS/images_for_markdown/blob/master/2016-11-29_CD_profile.jpg?raw=true)

中心差分的插值格式为

$$
\phi_e = \phi_C + \frac{\phi_E-\phi_C}{x_E-x_C}(x_e-x_C)
$$

也可写成

$$
\phi_e = \frac{x_E-x_e}{x_E-x_C} \phi_C + \frac{x_e-x_C}{x_E-x_C} \phi_E
$$

插值权函数`w(x)`为

$$
w(x) = \frac{x_E-x_e}{x_E-x_C} = \frac{d_{Ce}}{d_{Ce}+d_{eE}}
$$

三维情况下，权函数距离在面法向量上的投影之比。

`OpenFOAM`中中心差分格式的记法是`linear`，该类继承自 `surfaceInterpolationScheme<Type>`，重写了`weights()`和`linearInterpolate()`两个函数。`weights()`函数

```cpp
tmp<surfaceScalarField> weights
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return this->mesh().surfaceInterpolation::weights();
}
```
`weights()`函数调用了`surfaceInterpolation`类的`weights()`函数，该函数在[surfaceInterpolation.C](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolation.C)文件中，实现如下

```cpp
const Foam::surfaceScalarField&
Foam::surfaceInterpolation::weights() const
{
    if (!weights_)
    {
        makeWeights();
    }

    return (* weights_);
}

```

其中调用了`makeWeights()`函数

```cpp
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
}
```

### 3rd order QUICK scheme

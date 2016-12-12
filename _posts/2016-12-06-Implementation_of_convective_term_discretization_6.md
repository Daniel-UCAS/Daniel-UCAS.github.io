---
layout: post
title: Implementation of convective term discretization (6)
author: Daniel-UCAS
date: 2016-12-06
category: OpenFOAM
tags:
- OpenFOAM
- Convetive term
- Discretiazation
excerpt: Examples of limited schemes in OpenFOAM.
---

* TOC
{:toc}

### Examples of TVD schemes in OpenFOAM
上一篇文章说明了`OpenFOAM`中限制器函数的实现形式为
$$
(\phi_f)_{TVD} = \phi_P + \psi(r)(1-\omega_{CD})(\phi_N - \phi_P)
$$
下面看一下各种不同的限制器格式在`OpenFOAM`中的实现形式。
#### 1st order upwind scheme
##### 离散方法
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
##### 程序实现
继承自`limitedSurfaceInterpolationScheme<Type>`类，主要重写了`limiter`函数和两个`weights`函数

```cpp
virtual tmp<surfaceScalarField> limiter
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "upwindLimiter",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("upwindLimiter", dimless, 0.0)
        )
    );
}
tmp<surfaceScalarField> weights() const
{
    return pos(this->faceFlux_);
}
virtual tmp<surfaceScalarField> weights
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return weights();
}
```
### Examples of NVD schemes in OpenFOAM

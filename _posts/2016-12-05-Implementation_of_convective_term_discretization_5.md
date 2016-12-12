---
layout: post
title: Implementation of convective term discretization (5)
author: Daniel-UCAS
date: 2016-12-05
category: OpenFOAM
tags:
- OpenFOAM
- Convetive term
- Discretiazation
excerpt: High resolution schemes implementation in OpenFOAM.
---

* TOC
{:toc}

### High resolution scheme in unstructureed grid system

高阶/高分辨率格式通常需要三个及以上的节点作为模板，非结构网格的存储方式导致无法找到一个确定的远上游节点(U)，因此在程序编写时需要虚拟出一个U节点(如图所示)。
![virtual upwind node](https://github.com/Daniel-UCAS/images_for_markdown/blob/master/2016-12-01_virtual_upwind_node.jpg?raw=true)
在OpenFOAM中将上游节点C和下游节点D的连线反向延长，假定C位于$\overline{DU}$线段的中点，因而有

$$
\phi_D - \phi_U = \nabla \phi_C \cdot d_{UD} = 2\nabla \phi_C \cdot d_{CD}
$$
在上述假定下，可以推算出U点的变量值$\phi_U$.
$$
\phi_U = \phi_D -2\nabla \phi_C \cdot d_{CD}
$$
从而可以计算出TVD框架下的归一化变量r
$$
r = \frac{\phi_C-\phi_U}{\phi_D-\phi_C} = 2\frac{\nabla \phi_C \cdot d_{CD}}{\phi_D-\phi_C} - 1
$$
具体程序实现在[NVDTVD.H](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/NVDTVD.H)文件中的`r()`函数
```cpp
scalar r
(
    const scalar faceFlux,
    const scalar phiP,
    const scalar phiN,
    const vector& gradcP,
    const vector& gradcN,
    const vector& d
) const
{
    scalar gradf = phiN - phiP;

    scalar gradcf;

    if (faceFlux > 0)
    {
        gradcf = d & gradcP;
    }
    else
    {
        gradcf = d & gradcN;
    }

    if (mag(gradcf) >= 1000*mag(gradf))
    {
        return 2*1000*sign(gradcf)* sign(gradf) - 1;
    }
    else
    {
        return 2*(gradcf/gradf) - 1;
    }
}
};
```
TVD格式实现时，都会返回一个权重函数`weights()`，位于[`limitedSurfaceInterpolationScheme.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.C)。

[`Gamma`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/Gamma/Gamma.H)格式和[`SFCD`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/SFCD/SFCD.H)基于`phict()`函数计算，其定义为
$$
\tilde{\phi_C} = \frac{\phi_C-\phi_U}{\phi_D-\phi_U} = 1 - \frac{\phi_D-\phi_C}{2(\nabla\phi)_ C \cdot d}
$$
这一公式见于H.JASAK和H.G.WELLER的论文[`High resolution NVD differencing scheme for arbitrarily unstructured meshes`](http://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0363(19990930)31:2%3C431::AID-FLD884%3E3.0.CO;2-T/abstract)，即方程(20)。对应的函数实现为
```cpp
scalar phict
(
    const scalar faceFlux,
    const scalar phiP,
    const scalar phiN,
    const vector& gradcP,
    const vector& gradcN,
    const vector& d
) const
{
    scalar gradf = phiN - phiP;

    scalar gradcf;

    if (faceFlux > 0)
    {
        gradcf = d & gradcP;
    }
    else
    {
        gradcf = d & gradcN;
    }

    if (mag(gradf) >= 1000*mag(gradcf))
    {
        return 1 - 0.5*1000*sign(gradcf)* sign(gradf);
    }
    else
    {
        return 1 - 0.5*gradf/gradcf;
    }
}
```

每种不同的TVD格式对应于一个不同的`limiter()`函数，例如[`SuperBee`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/SuperBee/SuperBee.H)
```cpp
scalar limiter
(
    const scalar cdWeight,
    const scalar faceFlux,
    const typename LimiterFunc::phiType& phiP,
    const typename LimiterFunc::phiType& phiN,
    const typename LimiterFunc::gradPhiType& gradcP,
    const typename LimiterFunc::gradPhiType& gradcN,
    const vector& d
) const
{
    scalar r = LimiterFunc::r
    (
        faceFlux, phiP, phiN, gradcP, gradcN, d
    );

    return max(max(min(2*r, 1), min(r, 2)), 0);
}
```
其对应的限制器为
$$
\psi(r)=max(0,min(1,2r),min(2,r))
$$

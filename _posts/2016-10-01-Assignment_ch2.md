---
layout: post
title: Assignment for math
tags: [CFD, wind energy,OpenFOAM,Fluent,program]
excerpt: math assignment
---
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>
## 分离变量法作业
程瑜
工程热物理研究所
### 作业1

#### 方程组转化

$$
(1)\qquad
\begin{cases}
    \frac{\partial u}{\partial t} = a^2 \frac{\partial^2 u}{\partial x^2},0 < x < l , t > 0 \\\\
    u(x,0) = \phi(x) \\\\
    u(0,t) = 0,\frac{\partial u(l,t)}{\partial x} + hu(l,t) = 0
\end{cases}
$$

令$V(x,t)=u(x,t)-\\phi(x)$，方程组化为

$$
(2)\qquad
\begin{cases}
    \frac{\partial V}{\partial t} = a^2 \frac{\partial^2 u}{\partial x^2} - a^2\phi''(x),0 < x < l , t > 0  \\\\
    V(x,0) = 0 \\\\
    V(0,t) = -\phi(0), \frac{\partial V(l,t)}{\partial x} + hV(l,t) = -\phi'(l) - h\phi(l)
\end{cases}
$$

将方程(2)按叠加原理拆成两个方程组(3)和(4)
$$
(3)\qquad
\begin{cases}
    \frac{\partial V}{\partial t} = a^2 \frac{\partial^2 u}{\partial x^2},0 < x < l , t > 0  \\\\
    V(x,0) = 0 \\\\
    V(0,t) = -\phi(0), \frac{\partial V(l,t)}{\partial x} + hV(l,t) = -\phi'(l) - h\phi(l)
\end{cases}
$$

aaaa

$$
(4)\qquad
\begin{cases}
    \frac{\partial V}{\partial t} = a^2 \frac{\partial^2 u}{\partial x^2} - a^2\phi''(x),0 < x < l , t > 0  \\\\
    V(x,0) = 0 \\\\
    V(0,t) = -\phi(0), \frac{\partial V(l,t)}{\partial x} + hV(l,t) = 0
\end{cases}
$$

#### 求解方程(3)
利用分离变量法,令\\(V(x,t) = X(x)*T(t)\\), 方程简化为\\(XT'=a^2X''T\\)
即
$$
\frac{X''}{X} = a^2\frac{T}{T} = -\lambda
$$
方程(3)转化为
$$
\begin{cases}
T' + a^2\lambda T = 0 \\\\
T(0) = 0
\end{cases}
$$

##### 求解X(x)
$$x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}$$

#### 求解方程(4)

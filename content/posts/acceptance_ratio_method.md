---
title: "Acceptance_ratio_method"
date: 2021-12-09T19:47:59+08:00
# weight: 1
# tags: ["first"]
author: "Haomiao"
showToc: true
TocOpen: false
draft: false
hidemeta: false
comments: false
description: "BAR 方法推导"
disableHLJS: true # to disable highlightjs
disableShare: false
disableHLJS: false
hideSummary: false
searchHidden: true
ShowReadingTime: true
ShowBreadCrumbs: true
ShowPostNavLinks: true
math: true
# cover:
#     image: "<image path/url>" # image path/url
#     alt: "<alt text>" # alt text
#     caption: "<text>" # display caption under cover
#     relative: false # when using page bundles set this to true
#     hidden: true # only hide on current single page
# editPost:
#     URL: "https://github.com/<path_to_repo>/content"
#     Text: "Suggest Changes" # edit text
#     appendFilePath: true # to append file path to Edit link
---
## 问题定义

计算两个系统(0, 1) 在两个模拟中的自由能之差 $\Delta F$，在系统0, 1中的势能分别为 $U_0, U_1$.

1. 在这两个系统中的配分函数 分别为$Q_0,Q_1$  有如下关系
    $$
    \frac{Q_0}{Q_1} = \frac{Q_0}{Q_1} \frac{\int d\mathbf{r}^N w(\mathbf{r}^N) \exp(-\beta(U_0+U_1)))}{\int d\mathbf{r}^N w(\mathbf{r}^N) \exp(-\beta(U_0+U_1)))}
            = \frac{\left<w\exp(-\beta U_0)\right>_0}{\left<w\exp(-\beta U_1)\right>_1}
    $$

2. 自由能只差可以表示为配分函数的形式
    $$
    \beta\Delta F = \ln(Q_1/ Q_0)  = \ln\left<w\exp(-\beta U_0)\right>_1 - \ln \left<w\exp(-\beta U_1)\right>_0
    $$

3. 为了得到一个较好的 $\Delta F$ 的估计值，需要 $\Delta F$ 的方差最小，因此我们对 $\Delta F$的方差 (variance) 做一个估计
    $$
    \begin{aligned} \sigma_{\beta \Delta F}^{2}=& \frac{\left\langle\left[w \exp \left(-\beta \mathcal{U}_{1}\right)\right]^{2}\right\rangle_{0}-\left\langle w \exp \left(-\beta \mathcal{U}_{1}\right)\right\rangle_{0}^{2}}{n_{0}\left\langle w \exp \left(-\beta \mathcal{U}_{1}\right)\right\rangle_{0}^{2}} + \frac{\left\langle\left[w \exp \left(-\beta \mathcal{U}_{0}\right)\right]^{2}\right\rangle_{1}-\left\langle w \exp \left(-\beta \mathcal{U}_{0}\right)\right\rangle_{1}^{2}}{n_{1}\left\langle w \exp \left(-\beta \mathcal{U}_{0}\right)\right\rangle_{1}^{2}} \\
    =& \int \mathrm{d} \mathbf{r}^{N}\left\{\left[\left(Q_{0} / n_{0}\right) \exp \left(-\beta \mathcal{U}_{1}\right)+\left(Q_{1} / n_{1}\right) \exp \left(-\beta \mathcal{U}_{0}\right)\right]\right. \left.\times w^{2} \exp \left[-\beta\left(\mathcal{U}_{0}+\mathcal{U}_{1}\right)\right]\right\} 
    \\ &\times \frac{1}{\left\{\int d r^{N} w \exp \left[-\beta\left(\mathcal{U}_{0}+\mathcal{U}_{1}\right)\right]\right\}^{2}}-\frac{1}{n_{0}}-\frac{1}{n_{1}} .
    \end{aligned}
    $$
    ​	在得到这一步时用了方差的一些性质

    * 由于两个模拟互相独立，因此， $Var(X+Y) = Var(X) + Var(Y)$
    * 对于随机变量$X$, 存在公式 $Var(X) = E(X^2) - E(X)^2$
    * 如果变量有相同的方差 $\sigma^2$, 则对于随机变量 $X$, 则有 $Var(\bar{X}) = \frac{1}{n} + \frac{n-1}{n}\rho$, 其中 $\rho$是平均关联系数 (average correlation)，在我们的模拟中，我们认为采样是相互独立的，$\rho$ 的值为 $0$.
    *  对于任意函数，其方差可以表是为 $Var[g(X)] \approx (g'(\mu(x)))^2 \sigma^2_x$, 因此对于 log 函数有 $Var(\log(a)) \approx \sigma_a^2 / \mu_a^2 $
    * [Expected value and variance of log(a)](https://stats.stackexchange.com/questions/57715/expected-value-and-variance-of-loga)

4. 从3中可以看到， 如果我们给 $w$ 乘上一个常数因子，其值也不会变。因此，可以选择这样一个 $w$, 使得
    $$
    \int \mathrm{d} \mathbf{r}^{N} w \exp \left[-\beta\left(\mathcal{U}_{0}+\mathcal{U}_{1}\right)\right]=\text { constant }
    $$

5. 对 3 以及 4 的约束求导，使用拉格朗日乘子法，可以得到
   $$
   w=\frac{\text { constant }}{\left(Q_{0} / n_{0}\right) \exp \left(-\beta \mathcal{U}_{1}\right)+\left(Q_{1} / n_{1}\right) \exp \left(-\beta \mathcal{U}_{0}\right)}
   $$

 6. 将 5 带入1， 可得
    $$
    \frac{Q_{0}}{Q_{1}}=\frac{\left\langle\left\{1+\exp \left[\beta\left(\mathcal{U}_{0}-\mathcal{U}_{1}+C\right)\right]\right\}^{-1}\right\rangle}{\left\langle\left\{1+\exp \left[\beta\left(\mathcal{U}_{1}-\mathcal{U}_{0}-C\right)\right]\right\}^{-1}\right\rangle} \exp (\beta C)
    $$

 7. 其中 $\exp(\beta C) = (Q_0n_1)(Q_1n_0)$, 进一步的定义 Fermi-Dirac 函数 $f(x) = 1 / (1 + \exp(\beta x))$， 则有
    $$
    \frac{Q_{0}}{Q_{1}}=\frac{\left\langle f\left(\mathcal{U}_{0}-\mathcal{U}_{1}+C\right)\right\rangle_{1}}{\left\langle f\left(\mathcal{U}_{1}-\mathcal{U}_{0}-C\right)\right\rangle_{0}} \exp (\beta C)
    $$

8. $$
   \begin{array}{l}\left\langle f\left(\mathcal{U}_{0}-\mathcal{U}_{1}+\mathrm{C}\right)\right\rangle_{1}=\frac{1}{\mathrm{n}_{1}} \sum_{\mathrm{m}} \mathrm{f}_{\mathrm{m}}\left(\mathcal{U}_{0}-\mathcal{U}_{1}+\mathrm{C}\right) \\ \left\langle f\left(\mathcal{U}_{1}-\mathcal{U}_{0}-\mathrm{C}\right)\right\rangle_{0}=\frac{1}{n_{0}} \sum_{\mathrm{m}^{\prime}} \mathrm{f}_{\mathrm{m}^{\prime}}\left(\mathcal{U}_{1}-\mathcal{U}_{0}-\mathrm{C}\right),\end{array}
   $$

9. $$
   \beta \Delta F=\ln \frac{\sum_{1} f\left(\mathcal{U}_{0}-\mathcal{U}_{1}+\mathrm{C}\right)}{\sum_{0} f\left(\mathcal{U}_{1}-\mathcal{U}_{0}-\mathrm{C}\right)}-\ln \left(n_{1} / n_{0}\right)+\beta C,
   $$

10. $$
    \beta \Delta F=-\ln \left(n_{1} / n_{0}\right)+\beta C .
    $$

11. $$
    \sum_{m} f\left(\mathcal{U}_{0}-\mathcal{U}_{1}+C\right)=\sum_{m^{\prime}} f\left(\mathcal{U}_{1}-\mathcal{U}_{0}-C\right)
    $$

    
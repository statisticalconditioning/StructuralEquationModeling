# Measurement Model Scaling and Identification
William Murrah  
  







## Conventional Drawing and Labeling of SEM Model



### SEM nomenclature conventions

| name   | upper case   | lower case   | usage                                                            |
|--------|----------------|------------|------------------------------------------------------------------|
| Lambda |$\Lambda$     | $\lambda$    | Loading of a manifest indicator onto a latent construct          |
| Psi    | $\Psi$       | $\psi$       | residual variance/covariance of contruct when endogenous         |
| Theta  | $\Theta$     | $\theta$     | residual variance/covariance of indicators                       |
| Sigma  | $\Sigma$     | $\sigma$     | $\Sigma$ is the model implied variance/covariance matrix;        |
|        |              |              | $\sigma$ is standard deviation, $\sigma^2$ variance of indicator.|
|        |              |              | $\sigma$ can also be covariance of indicator                     |


<!--html_preserve--><div id="htmlwidget-8897" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-8897">{"x":{"diagram":"\ndigraph Basic {\nnode [shape = circle]\nConstruct_1;\nConstruct_2;\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\nindicator_4; indicator_5; indicator_6;\n# Edges\nConstruct_1 -> indicator_1 [label = <&nbsp;&lambda;<sub>1,1</sub>>];\nConstruct_1 -> indicator_2 [label = <&nbsp;&lambda;<sub>2,1</sub>>];\nConstruct_1 -> indicator_3 [label = <&nbsp;&lambda;<sub>3,1</sub>>];\nConstruct_1:n -> Construct_1:n [dir=both, label = <&psi;<sub>1,1</sub>>];\n\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1,1</sub>>];\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2,2</sub>>];\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3,3</sub>>];\n\nConstruct_2 -> indicator_4 [label = <&nbsp;&lambda;<sub>1,2</sub>>];\nConstruct_2 -> indicator_5 [label = <&nbsp;&lambda;<sub>2,2</sub>>];\nConstruct_2 -> indicator_6 [label = <&nbsp;&lambda;<sub>3,2</sub>>];\nConstruct_2:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,2</sub>>];\n\nindicator_4:s -> indicator_4:s [dir=both, label = <&theta;<sub>4,4</sub>>];\nindicator_5:s -> indicator_5:s [dir=both, label = <&theta;<sub>5,5</sub>>];\nindicator_6:s -> indicator_6:s [dir=both, label = <&theta;<sub>6,6</sub>>];\n\nConstruct_1:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,1</sub>>];\n\n{rank = same; Construct_1; Construct_2}\n{rank = same; indicator_1; indicator_2; indicator_3; \nindicator_4; indicator_5; indicator_6;}\n} \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->
\[ \textbf{$\Sigma$} = \left[ \begin{array}{cccccc}
\sigma_{1,1}^2 & \sigma_{1,2} & \sigma_{1,3} & \sigma_{1,4} & \sigma_{1,5} & \sigma_{1,6} \\
\sigma_{2,1} & \sigma_{2,2}^2 & \sigma_{2,3} & \sigma_{2,4} & \sigma_{2,5} & \sigma_{2,6} \\
\sigma_{3,1} & \sigma_{3,2} & \sigma_{3,3}^2 & \sigma_{3,4} & \sigma_{3,5} & \sigma_{3,6} \\
\sigma_{4,1} & \sigma_{4,2} & \sigma_{4,3} & \sigma_{4,4}^2 &  \sigma_{4,5} & \sigma_{4,6} \\
\sigma_{5,1} & \sigma_{5,2} & \sigma_{5,3} & \sigma_{5,4} & \sigma_{5,5}^2 & \sigma_{5,6} \\
\sigma_{6,1} & \sigma_{6,2} & \sigma_{6,3} & \sigma_{6,4} & \sigma_{6,5} & \sigma_{6,6}^2 
\end{array} \right],\] 

\[ \textbf{$\Lambda$} = \left[ \begin{array}{cc}
\lambda_{1,1} & \lambda_{1,2} \\
\lambda_{2,1} & \lambda_{2,2} \\
\lambda_{3,1} & \lambda_{3,2}
\end{array} \right],\] 

\[ \textbf{$\Psi$} = \left[ \begin{array}{cc}
\psi_{1,1} & \psi_{1,2} \\
\psi_{2,1} & \psi_{2,2}
\end{array} \right],\] 

\[ \textbf{$\Lambda^\prime$} = \left[ \begin{array}{cc}
\lambda_{1,1} & \lambda_{2,1} & \lambda_{3,1} \\
\lambda_{1,2} & \lambda_{2,2} & \lambda_{3,2}
\end{array} \right],\] 

\[ \textbf{$\Theta$} = \left[ \begin{array}{cccccc}
\theta_{1,1} & 0 & 0 & 0 & 0 & 0 \\
0 & \theta_{2,2} & 0 & 0 & 0 & 0 \\
0 & 0 & \theta_{3,3} & 0 & 0 & 0 \\
0 & 0 & 0 & \theta_{4,4} &  0 & 0 \\
0 & 0 & 0 & 0 & \theta_{5,5} & 0 \\
0 & 0 & 0 & 0 & 0 & \theta_{6,6} 
\end{array} \right].\] 

### Fundamental SEM equation

$$
\Sigma = \Lambda \Psi \Lambda' + \Theta \tag{1}
$$

## Unidimensional Measurement Model

<!--html_preserve--><div id="htmlwidget-4260" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-4260">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\n# Edges\nPositive -> Glad [label = <&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n\t{rank = same; Positive;}\n\t{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->


```r
library(lavaan)

##Prepare data with sufficient statisitics##
mymeans<-matrix(c(3.06893, 2.92590, 3.11013), ncol=3,nrow=1)
mysd<-c(0.84194,0.88934,0.83470)
mat <- c(1.00000,
         0.55226, 1.00000,
         0.56256, 0.66307, 1.00000)
mycor <- getCov(mat, lower = TRUE)

##Transform correlation matrix to covariance matrix using information above##
mycov <- mysd %*% t(mysd)

rownames(mycor) <-c( "Glad", "Cheerful", "Happy")
colnames(mycor) <-c( "Glad", "Cheerful", "Happy")

rownames(mycov) <-c( "Glad", "Cheerful", "Happy")
colnames(mycov) <-c( "Glad", "Cheerful", "Happy")
mynob<-823
```





<head><style>
        table {
              font-family: serif;
              text-align: right;}
        th {
              padding: 1px 1px 5px 5px;
	        }
        td {
             padding: 1px 1px 5px 5px; }
      </style></head><table align="center" style="border-collapse: collapse; caption-side:top; font-size:11pt;"><caption style="text-align:center;">Descriptive Statistics</caption><tr>
<th style="border-left: 0px solid black;background-color: #FFFFFF;border-top: 2px solid gray;border-bottom: 1px solid gray;">&nbsp;</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Glad</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Cheerful</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-right:0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Happy</th>
</tr>
<tr>
<td  style="border-left: 0px solid black; ">Glad</td>
<td align="right" style="border-left: 0px solid black;">1.00</td>
<td align="right" style="border-left: 0px solid black;">0.55</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;">0.56</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Cheerful</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.55</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">1.00</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.66</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Happy</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.66</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">1.00</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Mean</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">3.07</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">2.93</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">3.11</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">SD</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.84</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.89</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.83</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Var</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.79</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.70</td>
</tr>
<tr>
<td colspan="4" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>


### Measurement Model: 3 indicators

 using correlations only (instead of variance/covariance matirx)


```
Found more than one class "Model" in cache; using the first, from namespace 'MatrixModels'
```

```
lavaan (0.5-20) converged normally after   9 iterations

  Number of observations                           823

  Estimator                                         ML
  Minimum Function Test Statistic                0.000
  Degrees of freedom                                 0

Model test baseline model:

  Minimum Function Test Statistic                0.000
  Degrees of freedom                                 0
  P-value                                           NA

User model versus baseline model:

  Comparative Fit Index (CFI)                    1.000
  Tucker-Lewis Index (TLI)                       1.000

Loglikelihood and Information Criteria:

  Loglikelihood user model (H0)              -1070.768
  Loglikelihood unrestricted model (H1)      -1070.768

  Number of free parameters                          1
  Akaike (AIC)                                2143.536
  Bayesian (BIC)                              2148.249
  Sample-size adjusted Bayesian (BIC)         2145.074

Root Mean Square Error of Approximation:

  RMSEA                                          0.000
  90 Percent Confidence Interval          0.000  0.000
  P-value RMSEA <= 0.05                          1.000

Standardized Root Mean Square Residual:

  SRMR                                           0.000

Parameter Estimates:

  Information                                 Expected
  Standard Errors                             Standard

Latent Variables:
                   Estimate  Std.Err  Z-value  P(>|z|)
  Positive =~                                         
    Cheerful          1.000                           

Variances:
                   Estimate  Std.Err  Z-value  P(>|z|)
    Positive          0.790    0.039   20.285    0.000
    Cheerful          0.000                           
```


```r
cat(file = 'topics/2_MeasurementModel/2b_ConfirmatoryFactorAnalysis/mplus/l.cheer.out')
```

<!--html_preserve--><div id="htmlwidget-8570" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-8570">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\n# Edges\nPositive -> Glad [label = <&nbsp;&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&nbsp;&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&nbsp;&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n\t{rank = same; Positive;}\n\t{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->



```
Warning in lav_samplestats_from_moments(sample.cov = sample.cov,
sample.mean = sample.mean, : lavaan WARNING: sample covariance can not be
inverted
```

```
lavaan (0.5-20) converged normally after  99 iterations

  Number of observations                           823

  Estimator                                         ML
  Minimum Function Test Statistic                0.000
  Degrees of freedom                                 0

Parameter Estimates:

  Information                                 Expected
  Standard Errors                             Standard

Latent Variables:
                   Estimate  Std.Err  Z-value  P(>|z|)
  Positive =~                                         
    Glad              0.841    0.021   40.570    0.000
    Cheerful          0.889    0.022   40.570    0.000
    Happy             0.834    0.021   40.570    0.000

Variances:
                   Estimate  Std.Err  Z-value  P(>|z|)
    Glad              0.000    0.000   13.060    0.000
    Cheerful          0.000    0.000   12.213    0.000
    Happy             0.000    0.000   13.191    0.000
    Positive          1.000                           
```





<!--html_preserve--><div id="htmlwidget-5896" style="width:576px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-5896">{"x":{"diagram":"\ndigraph Basic {\nnode [shape = circle]\nConstruct_1;\nConstruct_2;\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\nindicator_4; indicator_5; indicator_6;\n# Edges\nConstruct_1 -> indicator_1 [label = <&nbsp;&lambda;<sub>1,1</sub>>];\nConstruct_1 -> indicator_2 [label = <&nbsp;&lambda;<sub>2,1</sub>>];\nConstruct_1 -> indicator_3 [label = <&nbsp;&lambda;<sub>3,1</sub>>];\nConstruct_1:n -> Construct_1:n [dir=both, label = <&psi;<sub>1,1</sub>>];\n\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1,1</sub>>];\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2,2</sub>>];\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3,3</sub>>];\n\nConstruct_2 -> indicator_4 [label = <&nbsp;&lambda;<sub>1,2</sub>>];\nConstruct_2 -> indicator_5 [label = <&nbsp;&lambda;<sub>2,2</sub>>];\nConstruct_2 -> indicator_6 [label = <&nbsp;&lambda;<sub>3,2</sub>>];\nConstruct_2:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,2</sub>>];\n\nindicator_4:s -> indicator_4:s [dir=both, label = <&theta;<sub>4,4</sub>>];\nindicator_5:s -> indicator_5:s [dir=both, label = <&theta;<sub>5,5</sub>>];\nindicator_6:s -> indicator_6:s [dir=both, label = <&theta;<sub>6,6</sub>>];\n\nConstruct_1:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,1</sub>>];\n\n{rank = same; Construct_1; Construct_2}\n{rank = same; indicator_1; indicator_2; indicator_3; \nindicator_4; indicator_5; indicator_6;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->
\[ \textbf{$\Sigma$} = \left[ \begin{array}{cccccc}
\sigma_{1,1}^2 & \sigma_{1,2} & \sigma_{1,3} & \sigma_{1,4} & \sigma_{1,5} & \sigma_{1,6} \\
\sigma_{2,1} & \sigma_{2,2}^2 & \sigma_{2,3} & \sigma_{2,4} & \sigma_{2,5} & \sigma_{2,6} \\
\sigma_{3,1} & \sigma_{3,2} & \sigma_{3,3}^2 & \sigma_{3,4} & \sigma_{3,5} & \sigma_{3,6} \\
\sigma_{4,1} & \sigma_{4,2} & \sigma_{4,3} & \sigma_{4,4}^2 &  \sigma_{4,5} & \sigma_{4,6} \\
\sigma_{5,1} & \sigma_{5,2} & \sigma_{5,3} & \sigma_{5,4} & \sigma_{5,5}^2 & \sigma_{5,6} \\
\sigma_{6,1} & \sigma_{6,2} & \sigma_{6,3} & \sigma_{6,4} & \sigma_{6,5} & \sigma_{6,6}^2 
\end{array} \right],\] 

\[ \textbf{$\Lambda$} = \left[ \begin{array}{cc}
\lambda_{1,1} & \lambda_{1,2} \\
\lambda_{2,1} & \lambda_{2,2} \\
\lambda_{3,1} & \lambda_{3,2}
\end{array} \right],\] 

\[ \textbf{$\Psi$} = \left[ \begin{array}{cc}
\psi_{1,1} & \psi_{1,2} \\
\psi_{2,1} & \psi_{2,2}
\end{array} \right],\] 

\[ \textbf{$\Lambda^\prime$} = \left[ \begin{array}{cc}
\lambda_{1,1} & \lambda_{2,1} & \lambda_{3,1} \\
\lambda_{1,2} & \lambda_{2,2} & \lambda_{3,2}
\end{array} \right],\] 

\[ \textbf{$\Theta$} = \left[ \begin{array}{cccccc}
\theta_{1,1} & 0 & 0 & 0 & 0 & 0 \\
0 & \theta_{2,2} & 0 & 0 & 0 & 0 \\
0 & 0 & \theta_{3,3} & 0 & 0 & 0 \\
0 & 0 & 0 & \theta_{4,4} &  0 & 0 \\
0 & 0 & 0 & 0 & \theta_{5,5} & 0 \\
0 & 0 & 0 & 0 & 0 & \theta_{6,6} 
\end{array} \right].\] 

### Latent Cheer with one indicator

```
# Mplus file
l.cheer.inp
```
<!--html_preserve--><div id="htmlwidget-2198" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-2198">{"x":{"diagram":"\ndigraph Cheer {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nCheer;\n\n# Edges\n\nPositive -> Cheer [label = <&nbsp;&lambda;>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N];\nCheer:s -> Cheer:s [dir=both, label = <&theta;>];\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->


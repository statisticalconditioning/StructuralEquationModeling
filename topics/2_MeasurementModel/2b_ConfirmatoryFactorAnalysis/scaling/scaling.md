---
title: "Measurement Model Scaling and Identification"
author: "William Murrah"
date: ''
output:
  html_document:
    fig_height: 3
    fig_width: 5
  pdf_document:
    fig_height: 3
    fig_width: 5
  word_document:
    fig_height: 3
    fig_width: 5
---







## Conventional Drawing and Labeling of SEM Model

<!--html_preserve--><div id="htmlwidget-5679" style="width:504px;height:504px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-5679">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nConstruct_1;\n\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\n\n# Edges\nConstruct_1 -> indicator_1 [label = <&nbsp;&lambda;<sub>1</sub>>];\nConstruct_1 -> indicator_2 [label = <&nbsp;&lambda;<sub>2</sub>>];\nConstruct_1 -> indicator_3 [label = <&nbsp;&lambda;<sub>3</sub>>];\nConstruct_1:n -> Construct_1:n [dir=both, label = <&psi;>]\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1</sub>>]\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2</sub>>]\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n\t{rank = same; Construct_1;}\n\t{rank = same; indicator_1; indicator_2; indicator_3;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->
## 

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
myvarcov <- outer(mysd, mysd, FUN="*")
mycov <- mycor * myvarcov

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

### Fundamental SEM equation
$$
\Sigma = \Lambda \Psi \Lambda + \Theta \tag{1}
$$


### Latent Cheer with one indicator



```
# Mplus file
l.cheer.inp
```
<!--html_preserve--><div id="htmlwidget-961" style="width:504px;height:504px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-961">{"x":{"diagram":"\ndigraph Cheer {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nCheer;\n\n# Edges\n\nPositive -> Cheer [label = <&nbsp;&lambda;>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N]\nCheer:s -> Cheer:s [dir=both, label = <&theta;>]\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->


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

<!--html_preserve--><div id="htmlwidget-5598" style="width:504px;height:504px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-5598">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\n# Edges\nPositive -> Glad [label = <&nbsp;&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&nbsp;&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&nbsp;&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n\t{rank = same; Positive;}\n\t{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->

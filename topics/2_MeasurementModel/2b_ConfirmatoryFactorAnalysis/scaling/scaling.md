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


<!--html_preserve--><div id="htmlwidget-3144" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-3144">{"x":{"diagram":"\ndigraph Basic {\nnode [shape = circle]\nConstruct;\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\n# Edges\nConstruct -> indicator_1 [label = <&nbsp;&lambda;<sub>1,1</sub>>];\nConstruct -> indicator_2 [label = <&nbsp;&lambda;<sub>2,1</sub>>];\nConstruct -> indicator_3 [label = <&nbsp;&lambda;<sub>3,1</sub>>];\nConstruct:n -> Construct:n [dir=both, label = <&psi;<sub>1,1</sub>>];\n\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1,1</sub>>];\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2,2</sub>>];\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3,3</sub>>];\n\n{rank = same; indicator_1; indicator_2; indicator_3;}\n} \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->
\[ \textbf{$\Sigma$} = \left[ \begin{array}{cccccc}
\sigma_{1,1}^2 & \sigma_{1,2} & \sigma_{1,3} \\
\sigma_{2,1} & \sigma_{2,2}^2 & \sigma_{2,3} \\
\sigma_{3,1} & \sigma_{3,2} & \sigma_{3,3}^2 
\end{array} \right],\] 

\[ \textbf{$\Lambda$} = \left[ \begin{array}{cc}
\lambda_{1,1} \\
\lambda_{2,1} \\
\lambda_{3,1}
\end{array} \right],\] 

\[ \textbf{$\Psi$} = \left[ \begin{array}{cc}
\psi_{1,1} 
\end{array} \right],\] 

\[ \textbf{$\Lambda^\prime$} = \left[ \begin{array}{cc}
\lambda_{1,1} & \lambda_{2,1} & \lambda_{3,1} 
\end{array} \right],\] 

\[ \textbf{$\Theta$} = \left[ \begin{array}{cccccc}
\theta_{1,1} & 0 & 0  \\
0 & \theta_{2,2} & 0  \\
0 & 0 & \theta_{3,3} 
\end{array} \right].\] 

### Fundamental SEM equation

$$
\Sigma = \Lambda \Psi \Lambda' + \Theta \tag{1}
$$

## Measurement Model: 3 indicators


<!--html_preserve--><div id="htmlwidget-233" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-233">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\n# Edges\nPositive -> Glad [label = <&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n\t{rank = same; Positive;}\n\t{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->


```r
library(lavaan)

##Prepare data with sufficient statisitics##
mymeans<-matrix(c(3.06893, 2.92590, 3.11013), ncol=3,nrow=1)
mysd<-c(0.84194,0.88934,0.83470)
mat <- c(1.00000,
         0.55226, 1.00000,
         0.56256, 0.60307, 1.00000)
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
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.60</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Happy</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.60</td>
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

### Mplus input file (`marker1.inp`)


```
 [1] TITLE: Positive Affect with Marker Variable Scaling    
 [2] DATA:                                                  
 [3]     FILE IS PAcorr2t.dat;                              
 [4]     TYPE IS CORRELATION;                               
 [5]     NOBSERVATIONS ARE 823;                             
 [6] VARIABLE:                                              
 [7]      NAMES ARE glad1 cheer1 happy1 glad2 cheer2 happy2;
 [8]      usevariables are glad1-happy1;                    
 [9] ANALYSIS:                                              
[10]                                                        
[11] MODEL:                                                 
[12]     Positive BY glad1 cheer1 happy1;                   
[13] OUTPUT: stdyx;                                         
```


### Mplus output using Correlation Table and first indicator fixed to 1.00


```
 [1] MODEL RESULTS                                                 
 [2]                                                     Two-Tailed
 [3]                     Estimate       S.E.  Est./S.E.    P-Value 
 [4]  POSITIVE BY                                                  
 [5]     GLAD1              1.000      0.000    999.000    999.000 
 [6]     CHEER1             1.072      0.061     17.612      0.000 
 [7]     HAPPY1             1.092      0.062     17.624      0.000 
 [8]  Variances                                                    
 [9]     POSITIVE           0.515      0.049     10.521      0.000 
[10]  Residual Variances                                           
[11]     GLAD1              0.484      0.033     14.547      0.000 
[12]     CHEER1             0.408      0.033     12.211      0.000 
[13]     HAPPY1             0.385      0.034     11.484      0.000 
```



```r
lambda = matrix(c(1.00, 1.072, 1.092), nrow = 3)
lambda
```

```
      [,1]
[1,] 1.000
[2,] 1.072
[3,] 1.092
```

```r
psi    = matrix(.515)
psi
```

```
      [,1]
[1,] 0.515
```

```r
t(lambda)
```

```
     [,1]  [,2]  [,3]
[1,]    1 1.072 1.092
```

```r
theta  = diag(c(.484, .408, .385))
theta
```

```
      [,1]  [,2]  [,3]
[1,] 0.484 0.000 0.000
[2,] 0.000 0.408 0.000
[3,] 0.000 0.000 0.385
```

```r
sigma = lambda %*% psi %*% t(lambda) + theta
round(sigma,2)
```

```
     [,1] [,2] [,3]
[1,] 1.00 0.55 0.56
[2,] 0.55 1.00 0.60
[3,] 0.56 0.60 1.00
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
      </style></head><table align="center" style="border-collapse: collapse; caption-side:top; font-size:11pt;"><caption style="text-align:center;"></caption><tr>
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
<td  style="border-left: 0px solid black; border-top: hidden;">Cheerful</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.55</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">1.00</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.60</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Happy</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.60</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">1.00</td>
</tr>
<tr>
<td colspan="4" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>

## Measurement Model: 6 indicators, 2 latent constructs

<!--html_preserve--><div id="htmlwidget-9174" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-9174">{"x":{"diagram":"\ndigraph Basic {\nnode [shape = circle]\nConstruct_1;\nConstruct_2;\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\nindicator_4; indicator_5; indicator_6;\n# Edges\nConstruct_1 -> indicator_1 [label = <&nbsp;&lambda;<sub>1,1</sub>>];\nConstruct_1 -> indicator_2 [label = <&nbsp;&lambda;<sub>2,1</sub>>];\nConstruct_1 -> indicator_3 [label = <&nbsp;&lambda;<sub>3,1</sub>>];\nConstruct_1:n -> Construct_1:n [dir=both, label = <&psi;<sub>1,1</sub>>];\n\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1,1</sub>>];\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2,2</sub>>];\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3,3</sub>>];\n\nConstruct_2 -> indicator_4 [label = <&nbsp;&lambda;<sub>4,2</sub>>];\nConstruct_2 -> indicator_5 [label = <&nbsp;&lambda;<sub>5,2</sub>>];\nConstruct_2 -> indicator_6 [label = <&nbsp;&lambda;<sub>6,2</sub>>];\nConstruct_2:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,2</sub>>];\n\nindicator_4:s -> indicator_4:s [dir=both, label = <&theta;<sub>4,4</sub>>];\nindicator_5:s -> indicator_5:s [dir=both, label = <&theta;<sub>5,5</sub>>];\nindicator_6:s -> indicator_6:s [dir=both, label = <&theta;<sub>6,6</sub>>];\n\nConstruct_1:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,1</sub>>];\n\n{rank = same; Construct_1; Construct_2}\n{rank = same; indicator_1; indicator_2; indicator_3; \nindicator_4; indicator_5; indicator_6;}\n} \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->

\[ \textbf{$\Sigma$} = \left[ \begin{array}{cccccc}
\sigma_{1,1}^2 & \sigma_{1,2} & \sigma_{1,3} & \sigma_{1,4} & \sigma_{1,5} & \sigma_{1,6} \\
\sigma_{2,1} & \sigma_{2,2}^2 & \sigma_{2,3} & \sigma_{2,4} & \sigma_{2,5} & \sigma_{2,6} \\
\sigma_{3,1} & \sigma_{3,2} & \sigma_{3,3}^2 & \sigma_{3,4} & \sigma_{3,5} & \sigma_{3,6} \\
\sigma_{4,1} & \sigma_{4,2} & \sigma_{4,3} & \sigma_{4,4}^2 &  \sigma_{4,5} & \sigma_{4,6} \\
\sigma_{5,1} & \sigma_{5,2} & \sigma_{5,3} & \sigma_{5,4} & \sigma_{5,5}^2 & \sigma_{5,6} \\
\sigma_{6,1} & \sigma_{6,2} & \sigma_{6,3} & \sigma_{6,4} & \sigma_{6,5} & \sigma_{6,6}^2 
\end{array} \right],\] 

\[ \textbf{$\Lambda$} = \left[ \begin{array}{cc}
\lambda_{1,1} & 0 \\
\lambda_{2,1} & 0 \\
\lambda_{3,1} & 0 \\
0 & \lambda_{4,2} \\
0 & \lambda_{5,2} \\
0 & \lambda_{6,2}
\end{array} \right],\] 

\[ \textbf{$\Psi$} = \left[ \begin{array}{cc}
\psi_{1,1} & \psi_{1,2} \\
\psi_{2,1} & \psi_{2,2}
\end{array} \right],\] 

\[ \textbf{$\Lambda^\prime$} = \left[ \begin{array}{cc}
\lambda_{1,1} & \lambda_{2,1} & \lambda_{3,1} 
0 & 0 & 0\\
0 & 0 & 0 & \lambda_{4,2} & \lambda_{5,2} & \lambda_{6,2}
\end{array} \right],\] 

\[ \textbf{$\Theta$} = \left[ \begin{array}{cccccc}
\theta_{1,1} & 0 & 0 & 0 & 0 & 0 \\
0 & \theta_{2,2} & 0 & 0 & 0 & 0 \\
0 & 0 & \theta_{3,3} & 0 & 0 & 0 \\
0 & 0 & 0 & \theta_{4,4} &  0 & 0 \\
0 & 0 & 0 & 0 & \theta_{5,5} & 0 \\
0 & 0 & 0 & 0 & 0 & \theta_{6,6} 
\end{array} \right].\] 

```r
lambda = matrix(c(.712, .788, .768, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, .729, .764, .778), nrow = 6)
lambda
```

```
      [,1]  [,2]
[1,] 0.712 0.000
[2,] 0.788 0.000
[3,] 0.768 0.000
[4,] 0.000 0.729
[5,] 0.000 0.764
[6,] 0.000 0.778
```

```r
psi    = matrix(c(1.00, 0.561,
                  0.561, 1.00), nrow = 2)
psi
```

```
      [,1]  [,2]
[1,] 1.000 0.561
[2,] 0.561 1.000
```

```r
theta = diag(c(.491, .378, .409, .467, .416, .394))
theta
```

```
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 0.491 0.000 0.000 0.000 0.000 0.000
[2,] 0.000 0.378 0.000 0.000 0.000 0.000
[3,] 0.000 0.000 0.409 0.000 0.000 0.000
[4,] 0.000 0.000 0.000 0.467 0.000 0.000
[5,] 0.000 0.000 0.000 0.000 0.416 0.000
[6,] 0.000 0.000 0.000 0.000 0.000 0.394
```

```r
sigma = lambda %*% psi %*% t(lambda) + theta
round(sigma, 2)
```

```
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,] 1.00 0.56 0.55 0.29 0.31 0.31
[2,] 0.56 1.00 0.61 0.32 0.34 0.34
[3,] 0.55 0.61 1.00 0.31 0.33 0.34
[4,] 0.29 0.32 0.31 1.00 0.56 0.57
[5,] 0.31 0.34 0.33 0.56 1.00 0.59
[6,] 0.31 0.34 0.34 0.57 0.59 1.00
```


```r
vec <- dput(scan('topics/2_MeasurementModel/2b_ConfirmatoryFactorAnalysis/mplus/PAcorr2t.dat', quiet = TRUE))
```

c(1, 0.55226, 1, 0.56256, 0.60307, 1, 0.31889, 0.35898, 0.27757, 
1, 0.24363, 0.35798, 0.31889, 0.56014, 1, 0.32217, 0.36385, 0.32072, 
0.56164, 0.59738, 1)

```r
mat <- getCov(vec, lower = TRUE)
ztable(mat)
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
      </style></head><table align="center" style="border-collapse: collapse; caption-side:top; font-size:11pt;"><caption style="text-align:center;"></caption><tr>
<th style="border-left: 0px solid black;background-color: #FFFFFF;border-top: 2px solid gray;border-bottom: 1px solid gray;">&nbsp;</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">V1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">V2</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">V3</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">V4</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">V5</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-right:0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">V6</th>
</tr>
<tr>
<td  style="border-left: 0px solid black; ">V1</td>
<td align="right" style="border-left: 0px solid black;">1.00</td>
<td align="right" style="border-left: 0px solid black;">0.55</td>
<td align="right" style="border-left: 0px solid black;">0.56</td>
<td align="right" style="border-left: 0px solid black;">0.32</td>
<td align="right" style="border-left: 0px solid black;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;">0.32</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">V2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.55</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">1.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.60</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.36</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.36</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.36</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">V3</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.60</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">1.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.28</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.32</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.32</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">V4</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.32</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.36</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.28</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">1.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.56</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">V5</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.36</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.32</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">1.00</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.60</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">V6</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.32</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.36</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.32</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.56</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.60</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">1.00</td>
</tr>
<tr>
<td colspan="7" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>


## Adding Means to the Model

#### The means structure equation:
$$
\mu_y = \textbf{T} + \Lambda \textbf{A}. \tag{2}
$$


<!--html_preserve--><div id="htmlwidget-4712" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-4712">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\nnode [shape = triangle]\n\n1;\n# Edges\nPositive -> Glad [label = <&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n1-> Positive [label = 0];\n1 -> Glad [label = <&tau;<sub>1</sub>>];\n1 -> Cheer [label = <&tau;<sub>2</sub>>];\n1 -> Happy [label = <&tau;<sub>3</sub>>];\n\n{rank = same; Positive; 1;}\n{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->
















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


<!--html_preserve--><div id="htmlwidget-4102" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-4102">{"x":{"diagram":"\ndigraph Basic {\nnode [shape = circle]\nConstruct;\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\n# Edges\nConstruct -> indicator_1 [label = <&nbsp;&lambda;<sub>1,1</sub>>];\nConstruct -> indicator_2 [label = <&nbsp;&lambda;<sub>2,1</sub>>];\nConstruct -> indicator_3 [label = <&nbsp;&lambda;<sub>3,1</sub>>];\nConstruct:n -> Construct:n [dir=both, label = <&psi;<sub>1,1</sub>>];\n\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1,1</sub>>];\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2,2</sub>>];\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3,3</sub>>];\n\n{rank = same; indicator_1; indicator_2; indicator_3;}\n} \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->
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


<!--html_preserve--><div id="htmlwidget-1657" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-1657">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\n# Edges\nPositive -> Glad [label = <&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n\t{rank = same; Positive;}\n\t{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->


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

<!--html_preserve--><div id="htmlwidget-8569" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-8569">{"x":{"diagram":"\ndigraph Basic {\nnode [shape = circle]\nConstruct_1;\nConstruct_2;\nnode [shape = box]\nindicator_1; indicator_2; indicator_3;\nindicator_4; indicator_5; indicator_6;\n# Edges\nConstruct_1 -> indicator_1 [label = <&nbsp;&lambda;<sub>1,1</sub>>];\nConstruct_1 -> indicator_2 [label = <&nbsp;&lambda;<sub>2,1</sub>>];\nConstruct_1 -> indicator_3 [label = <&nbsp;&lambda;<sub>3,1</sub>>];\nConstruct_1:n -> Construct_1:n [dir=both, label = <&psi;<sub>1,1</sub>>];\n\nindicator_1:s -> indicator_1:s [dir=both, label = <&theta;<sub>1,1</sub>>];\nindicator_2:s -> indicator_2:s [dir=both, label = <&theta;<sub>2,2</sub>>];\nindicator_3:s -> indicator_3:s [dir=both, label = <&theta;<sub>3,3</sub>>];\n\nConstruct_2 -> indicator_4 [label = <&nbsp;&lambda;<sub>4,2</sub>>];\nConstruct_2 -> indicator_5 [label = <&nbsp;&lambda;<sub>5,2</sub>>];\nConstruct_2 -> indicator_6 [label = <&nbsp;&lambda;<sub>6,2</sub>>];\nConstruct_2:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,2</sub>>];\n\nindicator_4:s -> indicator_4:s [dir=both, label = <&theta;<sub>4,4</sub>>];\nindicator_5:s -> indicator_5:s [dir=both, label = <&theta;<sub>5,5</sub>>];\nindicator_6:s -> indicator_6:s [dir=both, label = <&theta;<sub>6,6</sub>>];\n\nConstruct_1:n -> Construct_2:n [dir=both, label = <&psi;<sub>2,1</sub>>];\n\n{rank = same; Construct_1; Construct_2}\n{rank = same; indicator_1; indicator_2; indicator_3; \nindicator_4; indicator_5; indicator_6;}\n} \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->

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


<!--html_preserve--><div id="htmlwidget-9887" style="width:480px;height:288px;" class="grViz"></div>
<script type="application/json" data-for="htmlwidget-9887">{"x":{"diagram":"\n    \ndigraph CFA {\n\nnode [shape = circle]\nPositive;\n\nnode [shape = box]\nGlad; Cheer; Happy;\n\nnode [shape = triangle]\n\n1;\n# Edges\nPositive -> Glad [label = <&lambda;<sub>1</sub>>];\nPositive -> Cheer [label = <&lambda;<sub>2</sub>>];\nPositive -> Happy [label = <&lambda;<sub>3</sub>>];\nPositive:n -> Positive:n [dir=both, label = <&psi;>,position = N]\nGlad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]\nCheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]\nHappy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]\n\n1-> Positive [label = <&alpha;> ];\n1 -> Glad [label = <&tau;<sub>1</sub>>];\n1 -> Cheer [label = <&tau;<sub>2</sub>>];\n1 -> Happy [label = <&tau;<sub>3</sub>>];\n\n{rank = same; Positive; 1;}\n{rank = same; Glad; Cheer; Happy;}\n}  \n","config":{"engine":"dot","options":null}},"evals":[]}</script><!--/html_preserve-->


### Mplus: `factorMeans.inp`

```
 [1] TITLE: Positive Affect with Marker Variable Scaling and Mean Structure
 [2] DATA:                                                                 
 [3]     FILE IS PAmeansdcorr2t.dat;                                       
 [4]     TYPE IS MEANS STDEVIATIONS CORRELATION;                           
 [5]     NOBSERVATIONS ARE 823;                                            
 [6] VARIABLE:                                                             
 [7]      NAMES ARE glad1 cheer1 happy1 glad2 cheer2 happy2;               
 [8]      usevariables are glad1-happy1;                                   
 [9] ANALYSIS:                                                             
[10]                                                                       
[11] MODEL:                                                                
[12]     PA BY glad1* cheer1 happy1;                                       
[13]     PA@1;                                                             
[14]     [PA@0];                                                           
```


### Mplus: `factorMeans.out`

```
 [1]                                                     Two-Tailed
 [2]                     Estimate       S.E.  Est./S.E.    P-Value 
 [3]  PA       BY                                                  
 [4]     GLAD1              0.604      0.029     21.041      0.000 
 [5]     CHEER1             0.684      0.030     22.631      0.000 
 [6]     HAPPY1             0.654      0.028     23.074      0.000 
 [7]  Means                                                        
 [8]     PA                 0.000      0.000    999.000    999.000 
 [9]  Intercepts                                                   
[10]     GLAD1              3.069      0.029    104.631      0.000 
[11]     CHEER1             2.926      0.031     94.438      0.000 
[12]     HAPPY1             3.110      0.029    106.954      0.000 
[13]  Variances                                                    
[14]     PA                 1.000      0.000    999.000    999.000 
[15]  Residual Variances                                           
[16]     GLAD1              0.343      0.024     14.550      0.000 
[17]     CHEER1             0.322      0.026     12.209      0.000 
[18]     HAPPY1             0.268      0.023     11.483      0.000 
```



## Longitudinal Measurement Model



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
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Glad1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Cheer1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Happy1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Glad2</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Cheer2</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-right:0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Happy2</th>
</tr>
<tr>
<td  style="border-left: 0px solid black; ">Glad1</td>
<td align="right" style="border-left: 0px solid black;">0.71</td>
<td align="right" style="border-left: 0px solid black;">0.41</td>
<td align="right" style="border-left: 0px solid black;">0.40</td>
<td align="right" style="border-left: 0px solid black;">0.23</td>
<td align="right" style="border-left: 0px solid black;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;">0.23</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Cheer1</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.41</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.79</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.45</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.29</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Happy1</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.40</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.45</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.70</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.22</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Glad2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.23</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.43</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.40</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Cheer2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.29</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.43</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.83</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.46</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Happy2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.23</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.22</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.40</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.46</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.71</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Mean</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">3.07</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">2.93</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">3.11</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">3.03</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">2.86</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">3.09</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">SD</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.84</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.89</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.83</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.84</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.91</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.84</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Var</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.79</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.70</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.83</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.70</td>
</tr>
<tr>
<td colspan="7" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>


### Mplus: `effectCodingMeansLongitudinal.inp`

```
 [1] TITLE:                                                       
 [2]     Effects coding with means and Longitudinal               
 [3]     with 2 timepoints                                        
 [4] DATA:                                                        
 [5]     FILE IS PAmeansdcorr2t.dat;                              
 [6]     TYPE IS MEANS STD CORRELATION;                           
 [7]     NOBSERVATIONS ARE 823;                                   
 [8] VARIABLE:                                                    
 [9]     NAMES ARE Glad1 Cheer1 Happy1 Glad2 Cheer2 Happy2;       
[10]     !Note: glad1 cheer1 and happy1 are from Time 1           
[11]     !glad2 cheer2 and happy2 are from time 2                 
[12]                                                              
[13]     USEVARIABLES ARE Glad1 Cheer1 Happy1 Glad2 Cheer2 Happy2;
[14] MODEL:                                                       
[15]     Pos1 BY Glad1*(L1)                                       
[16]            Cheer1(L2)                                        
[17]            Happy1(L3) ; !Label for constraints               
[18]     Pos2 BY Glad2*(L4)                                       
[19]            Cheer2(L5)                                        
[20]            Happy2(L6);                                       
[21]     [Pos1]; !estimate latent means                           
[22]     [Pos2];                                                  
[23]     Glad1 with Glad2; !allow correlated residuals across time
[24]     Cheer1 with Cheer2;                                      
[25]     Happy1 with Happy2;                                      
[26]     [Glad1](T1)                                              
[27]     [Cheer1](T2)                                             
[28]     [Happy1](T3); !estimate and constrain the intercepts     
[29]     [Glad2](T4)                                              
[30]     [Cheer2](T5)                                             
[31]     [Happy2](T6);                                            
[32]     MODEL CONSTRAINT: L1= 3 - L2 - L3;                       
[33]                       L4 = 3 - L5 - L6;                      
[34]                       T1= 0 - T2 - T3;                       
[35]                       T4 = 0 - T5 - T6;                      
[36] OUTPUT:                                                      
[37]     TECH1;                                                   
[38]     RESIDUAL;                                                
[39]     STANDARDIZED;                                            
```

### Mplus: `effectcodingmeanslongitudinal.out`

```
 [1] MODEL RESULTS                                                 
 [2]                                                     Two-Tailed
 [3]                     Estimate       S.E.  Est./S.E.    P-Value 
 [4]  POS1     BY                                                  
 [5]     GLAD1              0.925      0.031     29.969      0.000 
 [6]     CHEER1             1.075      0.032     33.613      0.000 
 [7]     HAPPY1             0.999      0.031     32.008      0.000 
 [8]  POS2     BY                                                  
 [9]     GLAD2              0.934      0.031     30.394      0.000 
[10]     CHEER2             1.061      0.032     33.062      0.000 
[11]     HAPPY2             1.004      0.031     32.116      0.000 
[12]  POS2     WITH                                                
[13]     POS1               0.234      0.020     11.509      0.000 
[14]  GLAD1    WITH                                                
[15]     GLAD2              0.031      0.015      2.055      0.040 
[16]  CHEER1   WITH                                                
[17]     CHEER2             0.018      0.016      1.076      0.282 
[18]  HAPPY1   WITH                                                
[19]     HAPPY2            -0.011      0.014     -0.767      0.443 
[20]  Means                                                        
[21]     POS1               3.035      0.025    120.379      0.000 
[22]     POS2               2.992      0.025    117.661      0.000 
[23]  Intercepts                                                   
[24]     GLAD1              0.260      0.095      2.737      0.006 
[25]     CHEER1            -0.338      0.098     -3.435      0.001 
[26]     HAPPY1             0.078      0.096      0.808      0.419 
[27]     GLAD2              0.231      0.093      2.469      0.014 
[28]     CHEER2            -0.319      0.097     -3.275      0.001 
[29]     HAPPY2             0.088      0.095      0.933      0.351 
[30]  Variances                                                    
[31]     POS1               0.420      0.026     16.126      0.000 
[32]     POS2               0.426      0.026     16.106      0.000 
[33]  Residual Variances                                           
[34]     GLAD1              0.349      0.023     15.214      0.000 
[35]     CHEER1             0.306      0.025     12.199      0.000 
[36]     HAPPY1             0.278      0.022     12.455      0.000 
[37]     GLAD2              0.335      0.023     14.746      0.000 
[38]     CHEER2             0.343      0.026     13.084      0.000 
[39]     HAPPY2             0.274      0.023     12.153      0.000 
```



```r
lambda <- matrix(c(.925, 1.075, .999, 0, 0, 0,
                   0, 0, 0, .934, 1.061, 1.004), nrow = 6)
lambda
```

```
      [,1]  [,2]
[1,] 0.925 0.000
[2,] 1.075 0.000
[3,] 0.999 0.000
[4,] 0.000 0.934
[5,] 0.000 1.061
[6,] 0.000 1.004
```


```r
psi <- matrix(c(.420, .234,
                .234, .426), nrow=2)
psi
```

```
      [,1]  [,2]
[1,] 0.420 0.234
[2,] 0.234 0.426
```


```r
t(lambda)
```

```
      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 0.925 1.075 0.999 0.000 0.000 0.000
[2,] 0.000 0.000 0.000 0.934 1.061 1.004
```


```r
theta <- diag(c(.349, .306, .278, .335, .343, .274))
theta[1 ,4] <- theta[4, 1] <- .031 
theta[2, 5] <- theta[5, 2] <- .018
theta[3, 6] <- theta[6, 3] <- -.011
theta
```

```
      [,1]  [,2]   [,3]  [,4]  [,5]   [,6]
[1,] 0.349 0.000  0.000 0.031 0.000  0.000
[2,] 0.000 0.306  0.000 0.000 0.018  0.000
[3,] 0.000 0.000  0.278 0.000 0.000 -0.011
[4,] 0.031 0.000  0.000 0.335 0.000  0.000
[5,] 0.000 0.018  0.000 0.000 0.343  0.000
[6,] 0.000 0.000 -0.011 0.000 0.000  0.274
```



```r
tau <- matrix(c(.260, -.338, .078, .231, -.319, .088),nrow = 6)
tau
```

```
       [,1]
[1,]  0.260
[2,] -0.338
[3,]  0.078
[4,]  0.231
[5,] -0.319
[6,]  0.088
```


```r
alpha <- matrix(c(3.035, 2.992), nrow=2)
alpha
```

```
      [,1]
[1,] 3.035
[2,] 2.992
```


```r
sigma <- lambda %*% psi %*% t(lambda) + theta
sigma
```

```
          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
[1,] 0.7083625 0.4176375 0.3881115 0.2331643 0.2296535 0.2173158
[2,] 0.4176375 0.7913625 0.4510485 0.2349477 0.2848946 0.2525562
[3,] 0.3881115 0.4510485 0.6971604 0.2183374 0.2480257 0.2237011
[4,] 0.2331643 0.2349477 0.2183374 0.7066237 0.4221549 0.3994755
[5,] 0.2296535 0.2848946 0.2480257 0.4221549 0.8225571 0.4537939
[6,] 0.2173158 0.2525562 0.2237011 0.3994755 0.4537939 0.7034148
```


```r
mu <- tau + lambda %*% alpha
mu
```

```
         [,1]
[1,] 3.067375
[2,] 2.924625
[3,] 3.109965
[4,] 3.025528
[5,] 2.855512
[6,] 3.091968
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
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Glad1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Cheer1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Happy1</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Glad2</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Cheer2</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-right:0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Happy2</th>
</tr>
<tr>
<td  style="border-left: 0px solid black; ">Glad1</td>
<td align="right" style="border-left: 0px solid black;">0.71</td>
<td align="right" style="border-left: 0px solid black;">0.41</td>
<td align="right" style="border-left: 0px solid black;">0.40</td>
<td align="right" style="border-left: 0px solid black;">0.23</td>
<td align="right" style="border-left: 0px solid black;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;">0.23</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Cheer1</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.41</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.79</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.45</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.29</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Happy1</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.40</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.45</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.70</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.22</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Glad2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.23</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.43</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.40</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Cheer2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.29</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.43</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.83</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.46</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">Happy2</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.23</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.27</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.22</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.40</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.46</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.71</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Mean</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">3.07</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">2.93</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">3.11</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">3.03</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">2.86</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">3.09</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;background-color: #E5E4E2; ">SD</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.84</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.89</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.83</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.84</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #E5E4E2;">0.91</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E5E4E2;">0.84</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Var</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.79</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.70</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.71</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;">0.83</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;">0.70</td>
</tr>
<tr>
<td colspan="7" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>


















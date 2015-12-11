library(DiagrammeR)

dt <- grViz("
    
digraph CFA {

node [shape = circle]
Positive;

node [shape = box]
Glad; Cheer; Happy;

# Edges
Positive -> Glad [label = <&lambda;<sub>1</sub>>];
Positive -> Cheer [label = <&lambda;<sub>2</sub>>];
Positive -> Happy [label = <&lambda;<sub>3</sub>>];
Positive:n -> Positive:n [dir=both, label = <&psi;>,position = N]
Glad:s -> Glad:s [dir=both, label = <&theta;<sub>1</sub>>]
Cheer:s -> Cheer:s [dir=both, label = <&theta;<sub>2</sub>>]
Happy:s -> Happy:s [dir=both, label = <&theta;<sub>3</sub>>]

	{rank = same; Positive;}
	{rank = same; Glad; Cheer; Happy;}
}  
")

gr <- create_graph()
render_graph("topics/2_MeasurementModel/2b_ConfirmatoryFactorAnalysis/scaling/graph.gv")

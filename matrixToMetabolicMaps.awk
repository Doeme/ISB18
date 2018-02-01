NR>1{
	file="metabolic_maps/" $1 ".dot";
	printf("digraph{\nrankdir=LR;\n\"%s\"[shape=circle];\nnode[shape=plain];\n",$1) > file;
	i=2;
	while(length($i)>0){
		if($i<0){
			printf("v%d -> \"%s\";\n", i-1, $1)>file;
		}
		else if($i>0){
			printf("\"%s\" -> v%d;\n",$1,i-1)>file;
		}
		i++;
	}
	printf("}\n")>file;
}

BEGIN{
	print("digraph{");
}
NR==1{
	i=1;
	while(length($i)>0){
		printf("\"%s\";\n",$i);
		reactions[i]=$i;
		i++;
	}
	a=1;
}

function edge(from,to){
		printf("\"%s\" -> \"%s\";\n", from,to);
}

NR>1{
	metabolites[a]=$1;
	i=2;
	while(length($i)>0){
		if($i>0){
			edge(reactions[i-1],$1);
		}
		else if($i<0){
			edge($1,reactions[i-1]);
		}
		i++;
	}
	a++;
}
END{
	print("}");
}

NR==1{
	i=1;
	while(length($i)>0){
		reactions[i]=$i;
		i++;
	}
}
function sign(i){
	if(i>=0)
		return "+";
	else if(i<0)
		return "-";
}
function abs(i){
	if(i<0)
		return -i;
	return i;
}
NR>1{
	printf("\\frac{\\dd \\mathrm{%s}}{\\dd t} &= ", $1);
	i=2;
	first=(1==1);
	while(length($i)>0){
		if($i!=0){
			if(!first  || $i<0){
				printf(" %s", sign($i));
			}
			if(abs($i)>1){
				printf("%d ", abs($i));
			}
			printf("%s", reactions[i-1]);
			first=(1==0);

		}
		i++
	}
	printf("\\\\\n");
}

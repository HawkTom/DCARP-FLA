for ((j=2;j<=9;j++)){
for ((i=1; i <= 64; i++)){
	sshpass -p tongh1995 scp hao@10.16.28.227:/home/hao/fla/result/best/egl-e2-A-$j/instance$i.txt egl-e2-A-$j/instance$i.txt
	echo $i $j
}
}

#for ((i=2;i<=9;i++)){
#	mkdir egl-e2-A-$i
#}

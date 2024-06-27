ls ./raw/*_1.fq.gz | while read id
echo $id
do 
id2=${id%%_1.fq.gz}"_2.fq.gz" # '%%' 去除_1.fq.gz
trim_galore --q 25 -j 10 --stringency 4 --length 36 --paired $id $id2 --out ./trim # --q quality 默认（Default）20 ; -j thread ; 
#--stringency  指定修剪时接头序列的匹配严格度。接头序列的末端至少要有4个碱基与读取序列匹配，才会进行修剪。 Default 1 （可以适当放宽）
#--length 小于36的序列被去除

done

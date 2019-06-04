## Sources

目前最有用的Ref:

*  https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247487757&idx=1&sn=2fe6bc45acbd082a361b3d2c0d691563&chksm=fa46d430cd315d264a25b28ff45cb68c435f20800479381c29ab8d04b7d1b8eedc50731139b5&mpshare=1&scene=1&srcid=&key=bbe6b9114c05185534ac47902498bc54732f8cb939c47219da8d308e769b82fb8b4aa7f5514cadf4a639434eb753f5397c55ba597fe8e8ae45181bd722b605105648238cf88e929ec6a5398743ebc368&ascene=1&uin=NzA3NTE3MTMz&devicetype=Windows+10&version=62060739&lang=en&pass_ticket=SXN23ZaYL7dbyEsZV%2B30i3xdZWKqVM%2FnKEBdVugc9PNhHsL%2FzwllRcEbFZQerh4O

* https://mp.weixin.qq.com/s?__biz=MzU4NjU4ODQ2MQ==&mid=2247486039&idx=2&sn=887b456b6b762b75094eb70cfbad1ac6&chksm=fdf84215ca8fcb0308213948c32cf7a70449e2b5a0c2fcee4f6b80f650738ec1e00258934c96&mpshare=1&scene=1&srcid=&key=c8c9cb9453e0935031f4a12f7d3a5784c4026682a6cf81d6b04d42e00e93b2da51c4037221d45e851a0e91480113444778b6596784c8cbf3a64a28976ac19bc67a12a323c5af0283a5f16be6e2ec5038&ascene=1&uin=NzA3NTE3MTMz&devicetype=Windows+10&version=62060739&lang=en&pass_ticket=qsvmCxoFviqMD1c7jgwYAO3C7Uu8eH2USUAm%2F8ynoPIf5aNtlWgXD2Y1sNQLgw4I 

Ref:

*  https://github.com/stephenturner/oneliners
*  https://github.com/stephenturner/oneliners/blob/master/README.md#contents  

## Commands

### 涉及计算

[datamash](<https://www.gnu.org/software/datamash/>)

### 常用

输出第七列满足正则表达式的行：

```bash
awk '$7  ~ /^[a-f]/' file.txt
```

根据第二列输出`file.txt`文件中的唯一条目：

```
awk '!arr[$2]++' file.txt
```

> 以`awk '!arr[$0]++'`为例来拆解一下步骤：
>
> ```
> $ cat file.txt
> A23     B66     1234
> C56     D34     2334
> A23     B66     1234
> C56     D34     2334
> $ awk '{ print arr[$0]++}' file.txt
> 0
> 0
> 1
> 1
> $ awk '!arr[$0]++' file.txt
> A23     B66     1234
> C56     D34     2334
> ```
>
> `awk`首先抽取第一行，此时`arr[$0]`为`arr["A23 B66 1234"]`，值为空字符串，则有`!arr[$0]=1`，需要注意的是，`++`的值是在对`arr[$0]`值取反后才进行，最后打印本行。也就是说`!arr[$0]++`，其实是两个分开的并行动作，`!arr[$0]`是对关系表达式取反，`arr[$0]++`对`arr[$0]`进行`+1`的操作。
>
> 直到遇到相同的行，比如上面`file.txt`文件中第三行有重复的`A23 B66 1234`，在对第一行执行完操作后，`arr["A23 B66 1234"]++`后值为 1。对第二个`arr["A23 B66 1234"]++`行进行处理时，由于数组中已有此元素，则不再重新赋值，`arr["A23 B66 1234"]++`值为 1，对其取反，`!arr["A23 B66 1234"]=0`，关系表达式为假，则本行不再打印。

将`file.txt`文件中第一列的值相加：

```
awk '{sum+=$1} END {print sum}' file.txt
```

计算第二列的平均值：

```
awk '{x+=$2}END{print x/NR}' file.txt
```

将`file.txt`文件中所有的`foo`替换为 `bar`：

```
sed 's/foo/bar/g' file.txt
```

去除`file.txt`文件中开头的空格和制表符：

```
sed 's/^[ \t]*//' file.txt
```

去除`file.txt`文件中结尾的空格和制表符：

```
sed 's/[ \t]*$//' file.txt
```

去除`file.txt`文件中开头和结尾的空格和制表符：

```
sed 's/^[ \t]*//;s/[ \t]*$//' file.txt
```

删除`file.txt`文件中的空行：

```
sed '/^$/d' file.txt
```

删除`file.txt`文件中包括`EndOfUsefulData`及其之后的行：

```
sed -n '/EndOfUsefulData/,$!p' file.txt
#
# -n 只打印模式匹配的行
# /EndOfUsefulData/,$ 意思为从匹配 /EndOfUsefulData/ 的行到最后一行
# !p 打印不匹配的行
```

### Sequencing

fastq sequences length distribution  => 得到fq文件中序列长度的分布:  
```bash
zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  
```

reverse complement  => 反向互补:  

```
echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'
```

fastq2fasta:  

```
zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[1]\n";' | gzip -c > file.fasta.gz
```

split a multifasta file into single ones with csplit => fasta按>拆分:  

```
csplit -z -q -n 4 -f test test.fa /\>/ {*}
```
```
awk '/^>/{s=++d".fa"} {print > s}' multi.fa
```

single line fasta to multi-line of 50 characters in each line => 单行fa变多行:  

```
awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/50;i++) {for (j=1;j<=50;j++) printf "%s", $(i*50 +j); print ""}}' file
```

multi-line fasta to one-line => 一个多行fa文件变单行:  

```
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' file.fa
```
```
cat file.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}'
```

bam2bed:  

```
samtools view file.bam | perl -F'\t' -ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > file.bed
```
bam2wig:  

```
samtools mpileup -BQ0 file.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > file.wig.gz
```

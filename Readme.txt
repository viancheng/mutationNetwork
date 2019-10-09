突变网络分析模块用来根据maf文件分析群体中的基因协同和基因互斥事件；
用法（work.sh是一个示例）：
python mutationNetwork.py -i maf.list -k 3 -f 0.05 -n 30 -p1 0.05 -p2 0.01 -m vcfilter -fdr 1 -o outdir

参数说明及注意事项：
1.必须使用python3，建议/zfssz5/BC_PS/chengyuanfang/software/miniconda3/bin/python；
2.-i为输入文件，该文件的每一行是一个maf文件路径，请注意这些maf文件的header必须一致；必须设置；
3.-f为基因在群体中的最小频率，默认0.05，用户可自己设置，请注意若单个样本的频率超过f，则会报错；非必须设置；
4.-k为基因集大小，默认为2（即每行显示的是一对基因的协同/互斥结果），若需要一行显示多个基因的结果，则需要调整k值；
5.-n为maftools展示时的top n个基因，默认30，用户可自己设置；非必须设置；
6.-p1为协同事件的p-value显著性阈值，默认0.05，用户可自己设置；非必须设置；
7.-p2为互斥事件的p-value显著性阈值，默认0.01，用户可自己设置；非必须设置；
8.-m为是否过滤maf文件，vcfilter/vcall为两个可选项，前者会过滤maf仅保留常用突变类型，后者不进行过滤；非必须设置；
9.-fdr为false discover rate阈值（0-1之间），默认为1，值越大犯第一类错误可能性越大；非必须设置；
10.-o为输出目录，该目录不需要提前建立，输出结果为：outdir/cooccurrence.filtered.txt（协同事件）、outdir/exclusivity.filtered.txt（互
斥事件）、outdir/*.pdf（协同/互斥的图形展示），其他文件夹下为中间过程文件；必须设置；
11.如果需要更改paf图片，请修改outdir/maftools/runMaftools.r后重新运行；

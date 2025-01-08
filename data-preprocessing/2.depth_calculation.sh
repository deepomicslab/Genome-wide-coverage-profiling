f=$1; declare -i i=0;
b=$2
for((r=1;r<=22;r++));
#for r in $(samtools view $f | cut -f3 | uniq | sed "s/*//")
do
	samtools depth -J -s $f -r chr$r > $b.tmp.depth.$r.out &
	pids[$i]=$!
	regs[$i]=$r
	i=$i+1
done
for pid in ${pids[*]}; do wait $pid; done
for r in ${regs[*]}; do cat $b.tmp.depth.$r.out >> $b.depth.out; rm $b.tmp.depth.$r.out; done
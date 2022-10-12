for i in `seq 0 3`; do
  ssh nari@192.168.20.16$i ~/opttp/sync.sh
done

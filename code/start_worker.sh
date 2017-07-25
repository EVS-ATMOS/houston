echo "number of nodes"
export PBS_NUM_NODES=`wc -l ${PBS_NODEFILE} | cut -f1 -d" "`
echo $PBS_NUM_NODES
echo "===================="
echo "ssh into each node"


#====================================================================
scheduler_node=$1
hostIndex=0
echo "Pulling from schedule node:"
echo $scheduler_node
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
dask-worker --nprocs 1 --nthreads 16  $scheduler_node:8786

# # washington.n6c7
# echo ">>>================================<<<"
# echo ">>> benchmarking washington.n6c7"
# ./maxflow prf -f generate/washington.n6c7/washington.n6c7.max
# ./maxflow ppr -f generate/washington.n6c7/washington.n6c7.max -p 1
# ./maxflow ppr -f generate/washington.n6c7/washington.n6c7.max -p 2
# ./maxflow ppr -f generate/washington.n6c7/washington.n6c7.max -p 4
# ./maxflow ppr -f generate/washington.n6c7/washington.n6c7.max -p 8
# echo ""

# # liver.n6c10
# echo ">>>================================<<<"
# echo ">>> benchmarking liver.n6c10"
# ./maxflow prf -f generate/liver.n6c10/liver.n6c10.max
# ./maxflow ppr -f generate/liver.n6c10/liver.n6c10.max -p 1
# ./maxflow ppr -f generate/liver.n6c10/liver.n6c10.max -p 2
# ./maxflow ppr -f generate/liver.n6c10/liver.n6c10.max -p 4
# ./maxflow ppr -f generate/liver.n6c10/liver.n6c10.max -p 8
# echo ""

# # adhead.n6c10
# echo ">>>================================<<<"
# echo ">>> benchmarking adhead.n6c10"
# ./maxflow prf -f generate/adhead.n6c10/adhead.n6c10.max
# ./maxflow ppr -f generate/adhead.n6c10/adhead.n6c10.max -p 1
# ./maxflow ppr -f generate/adhead.n6c10/adhead.n6c10.max -p 2
# ./maxflow ppr -f generate/adhead.n6c10/adhead.n6c10.max -p 4
# ./maxflow ppr -f generate/adhead.n6c10/adhead.n6c10.max -p 8
# echo ""

# LB07-bunny-med
echo ">>>================================<<<"
echo ">>> benchmarking LB07-bunny-med"
./maxflow prf -f generate/LB07-bunny-med/LB07-bunny-med.max
./maxflow ppr -f generate/LB07-bunny-med/LB07-bunny-med.max -p 1
./maxflow ppr -f generate/LB07-bunny-med/LB07-bunny-med.max -p 2
./maxflow ppr -f generate/LB07-bunny-med/LB07-bunny-med.max -p 4
./maxflow ppr -f generate/LB07-bunny-med/LB07-bunny-med.max -p 8
echo ""
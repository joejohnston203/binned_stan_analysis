#!/bin/bash

outfile="parameters.txt"
#sims_path="../../../data/simulations_cuore/reduced_sims/combined_trees/"
sims_path="/nfs/cuore1/data/simulation/CUORE/2017/ntp/"

echo "parameters:" > $outfile

for file in `cat ../jags_inputs/ListMC.txt | sed 's/\(.*\)\.root.*/\1/'`; do
    echo ' - name: "p_'$file'"' >> $outfile
    echo '   lower_bound: 0.' >> $outfile
    echo '   upper_bound: 1000.' >> $outfile
    echo '   #prior: ""' >> $outfile
    echo '   shapes:' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["Ener2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==1&&MultipletIndex==1&&Layer==0"' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["Ener2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==1&&MultipletIndex==1&&Layer==1"' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["Ener2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==2&&MultipletIndex==1"' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["ESum2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==2&&MultipletIndex==1"' >> $outfile
done

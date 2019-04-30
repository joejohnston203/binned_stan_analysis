#!/bin/bash

# Note- A normal prior must be added manually for p_mu-R5m_cnaf127
#   prior: "normal(2.93e-4, 2.93e-5)"

outfile="parameters.txt"
#sims_path="../../../data/simulations_cuore/reduced_sims/combined_trees/"
sims_path="/nfs/cuore1/data/simulation/CUORE/2017/ntp/"

# Got parameter upper bound with:
# sed -n 's/upN[0-9]* <- \(.*\)/\1/p' Data4Jags.dat > upper_bounds.txt
upper_bounds_file="upper_bounds.txt"

# Combine file names and upper bounds
cat ../jags_inputs/ListMC.txt | sed 's/\(.*\)\.root.*/\1/' > temp_parameter_names.txt

echo "parameters:" > $outfile

paste temp_parameter_names.txt $upper_bounds_file | while read file upper_bound; do
    echo ' - name: "p_'$file'"' >> $outfile
    echo '   lower_bound: 0.' >> $outfile
    echo '   upper_bound: '$upper_bound >> $outfile
    echo '   #prior: ""' >> $outfile
    echo '   shapes:' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["Ener2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==1&&Layer==0"' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["Ener2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==1&&Layer==1"' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["Ener2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==2"' >> $outfile
    echo '    - path: "'$sims_path$file'.root"' >> $outfile
    echo '      renormalize: False' >> $outfile
    echo '      format: "root values"' >> $outfile
    echo '      tree: "outTree"' >> $outfile
    echo '      branches: ["ESum2"]' >> $outfile
    echo '      cut: "Included==1&&Multiplicity==2&&MultipletIndex==1"' >> $outfile
done

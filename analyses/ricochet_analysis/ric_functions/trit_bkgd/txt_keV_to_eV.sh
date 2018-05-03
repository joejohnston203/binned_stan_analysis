#!/bin/bash

sed -n 's/^\([0-9]*\)\.\([0-9]\{3\}\)\([0-9]*\)/\1\2.\3/p' ricochet_sig_back_recoils_keV.txt > ricochet_sig_back_recoils_eV.txt

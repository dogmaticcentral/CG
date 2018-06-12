
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL11.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 13+63+93+100+114+118+148+152+165
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 22+129
color red, red_clust
show spheres, red_clust


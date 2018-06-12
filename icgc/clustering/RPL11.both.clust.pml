
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL11.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 43+59+63+79+86+93+100+114+118+138
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 19+22+29+33+37+38+48+69+72+73+107+129+131
color red, red_clust
show spheres, red_clust

select blue_clust, resi 13+15+16+164+165+166+172
color blue, blue_clust
show spheres, blue_clust

select green_clust, resi 142+151+152
color green, green_clust
show spheres, green_clust

select cluster_3, resi 121+122
color orange, cluster_3
show spheres, cluster_3

deselect

select cluster_4, resi 125+127
color orange, cluster_4
show spheres, cluster_4

deselect

select cluster_5, resi 146+148
color orange, cluster_5
show spheres, cluster_5

deselect

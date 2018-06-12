
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL11.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 43+59+79+86+107+127+138+146+151
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 19+29+33+37+38+48+69+72+73+131
color red, red_clust
show spheres, red_clust

select blue_clust, resi 164+166+172
color blue, blue_clust
show spheres, blue_clust

select green_clust, resi 15+16
color green, green_clust
show spheres, green_clust

select cluster_3, resi 121+122
color orange, cluster_3
show spheres, cluster_3

deselect

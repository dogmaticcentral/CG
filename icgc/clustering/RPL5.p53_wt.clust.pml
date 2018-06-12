
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL5.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 3+78+82+137+174+258+287
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 91+94+97+201+208+226+236
color red, red_clust
show spheres, red_clust

select blue_clust, resi 107+169+170+248+251
color blue, blue_clust
show spheres, blue_clust

select green_clust, resi 184+185+189
color green, green_clust
show spheres, green_clust

select cluster_3, resi 54+147
color orange, cluster_3
show spheres, cluster_3

deselect

select cluster_4, resi 57+59
color orange, cluster_4
show spheres, cluster_4

deselect

select cluster_5, resi 68+71
color orange, cluster_5
show spheres, cluster_5

deselect

select cluster_6, resi 125+196
color orange, cluster_6
show spheres, cluster_6

deselect

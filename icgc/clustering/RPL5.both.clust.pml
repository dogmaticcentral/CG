
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL5.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 3+8+30+82+137+151+154+176+268+287+293
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 91+94+97+100+102+125+164+181+195+196+198+201+205+207+208+209+226+231+236
color red, red_clust
show spheres, red_clust

select blue_clust, resi 107+109+169+170+248+251
color blue, blue_clust
show spheres, blue_clust

select green_clust, resi 48+50+54+147
color green, green_clust
show spheres, green_clust

select cluster_3, resi 62+77+78
color orange, cluster_3
show spheres, cluster_3

deselect

select cluster_4, resi 184+185+189
color orange, cluster_4
show spheres, cluster_4

deselect

select cluster_5, resi 57+59
color orange, cluster_5
show spheres, cluster_5

deselect

select cluster_6, resi 68+71
color orange, cluster_6
show spheres, cluster_6

deselect

select cluster_7, resi 173+174
color orange, cluster_7
show spheres, cluster_7

deselect

select cluster_8, resi 257+258
color orange, cluster_8
show spheres, cluster_8

deselect

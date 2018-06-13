
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL5.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 78+82+83+87+99+102+103+104+106+107+108+110+127+153+156+160+164+165+166+168+169+170+179+180+182+195+247+248
color red, red_clust
show spheres, red_clust

select blue_clust, resi 22+23+26+27+28+29+33+35+36+39+47+48+50+53+54+63+71+72+146+147+149+150
color blue, blue_clust
show spheres, blue_clust

select green_clust, resi 91+92+94
color green, green_clust
show spheres, green_clust

select cluster_3, resi 207+208+223
color orange, cluster_3
show spheres, cluster_3

deselect

select cluster_4, resi 12+15
color orange, cluster_4
show spheres, cluster_4

deselect

select cluster_5, resi 174+175
color orange, cluster_5
show spheres, cluster_5

deselect

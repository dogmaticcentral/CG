
load /home/ivana/Dropbox/Sinisa/ribosomal/data/structures/5s-rRNP.RPL5.pdb

bg_color  white
hide all
show cartoon
color white

#isolated
select isolated, resi 8+30+62+109+151+154+173+176+231+257+268
show spheres, isolated
color gray, isolated

#clusters

select red_clust, resi 164+181+195+198
color red, red_clust
show spheres, red_clust

select blue_clust, resi 48+50
color blue, blue_clust
show spheres, blue_clust

select green_clust, resi 100+102
color green, green_clust
show spheres, green_clust


theoryCovness=0
simname=LightConeHector
grid=LightConeHectorPatchyWideHDvFFT_IR06

cd /exports/pierre/EFTofBOSS


for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.1
do
for number in 1
do
qsub -l nodes=node01:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.Prior$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.Prior$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 
done
done
done
done
done
done
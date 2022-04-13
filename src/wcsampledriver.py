
import rosypig as rp
from wc import WCField
from sampledriver import SampleDriver

class WCSampleDriver(SampleDriver):

    # New WC01 parameters
    n_time_steps     = 100
    use_drc_max_lz   = False

    # Irrelevant to WC01
    jacobi_max_iters = None
    sorDict          = None

    # Not implemented in MITgcm
    smooth2DDims     = None

    def range_approx_one(self):
        jid_list = []
        for Nx in self.NxList:
            for xi in self.xiList:

                _, write_dir, run_dir = self._get_dirs('range_approx_one',Nx,xi)

                self.smooth_writer(write_dir,
                                   xi=xi,
                                   smooth_apply=False,
                                   num_inputs=self.n_samples)

                matern = WCField(xdalike=self.mymodel,
                                 n_range=Nx,
                                 horizontal_factor=xi,
                                 isotropic=self.isotropic)
                matern.write_binaries(write_dir=write_dir,
                                      smoothOpNb=self.smoothOpNb)

                sim = rp.Simulation(name=f'{Nx:02}dx_{xi:02}xi_oi_ra1',
                                    run_dir=run_dir,
                                    obs_dir=write_dir,
                                    **self.dsim)

                # launch job
                sim.link_to_run_dir()

                sim.write_slurm_script()
                jid = sim.submit_slurm(**self.slurm)
                jid_list.append(jid)


    def smooth_writer(self, write_dir, xi, smooth_apply=True, num_inputs=1000):

        smooth = f"smooth{self.n_dims}D"
        alg = "WC01       "
        maskName = "mask"+self.gridloc

        file_contents = ' &SMOOTH_NML\n'+\
            f' smoothMdsDir = "smooth-output",\n'+\
            f' {smooth}Algorithm({self.smoothOpNb})=\'{alg}\',\n'+\
            f' {smooth}CreateOperator({self.smoothOpNb}) = .TRUE.,\n'+\
            f' {smooth}CalcNormFactor({self.smoothOpNb}) = .TRUE.,\n'+\
            f' {smooth}WriteSamples({self.smoothOpNb}) = .TRUE.,\n'+\
            f' {smooth}ConstHorizontal({self.smoothOpNb}) = .FALSE.,\n'+\
            f' {smooth}ConstVertical({self.smoothOpNb}) = .FALSE.,\n'+\
            f' {smooth}UseDRCMaxLz({self.smoothOpNb}) = .{str(self.use_drc_max_lz).upper()}.,\n'+\
            f' {smooth}Nbt({self.smoothOpNb}) = {self.n_time_steps},\n'+\
            f' {smooth}NbRand({self.smoothOpNb}) = {num_inputs},\n'+\
            f' {smooth}MaskName({self.smoothOpNb}) = "{maskName}",\n'+\
            ' &'
        fname = write_dir+f'/data.smooth'
        with open(fname,'w') as f:
            f.write(file_contents)


    def write_bash_script(self,stage,mysim):
        """Write a bash script for the next experiment stage
        """
        file_contents = '#!/bin/bash\n\n' +\
            f'#SBATCH -J {stage}\n' +\
            f'#SBATCH -o {stage}.%j.out\n' +\
            f'#SBATCH -e {stage}.%j.err\n' +\
            '#SBATCH -N 1\n' +\
            f'#SBATCH -n {mysim.procs_per_node}\n' +\
            f'#SBATCH -p {mysim.queue_name}\n' +\
            f'#SBATCH -t {mysim.time}\n'

        file_contents += f'\n\neval "$(conda shell.bash hook)"\n'+\
                f'conda activate {self.conda_env}\n\n'+\
                f'python3 -c '+\
                '"from wcsampledriver import WCSampleDriver;'+\
                f'oid = WCSampleDriver(\'{self.experiment}\',\'{stage}\')"\n'

        fname = self.dirs['main_run']+f'/submit_{self.experiment}.sh'
        with open(fname,'w') as f:
            f.write(file_contents)
        return fname

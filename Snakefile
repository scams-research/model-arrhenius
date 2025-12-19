for temp in [300, 350, 400, 500, 600, 700]:
    for length in [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]:
        rule:
            name:
                f"agcrse2_analysis_{temp}_{length}"
            input:
                'src/code/agcrse2_analysis.py',
                'src/code/mda_helper.py',
                [f'src/data/raw/agcrse2/agcrse2_T{temp}_{id}.lammpstrj' for id in range(1, 22)],
                'src/data/raw/agcrse2/structure.data'
            output:
                f'src/data/reduced/agcrse2/Dmsd{temp}_{length}.h5',
            conda:
                'environment.yml'
            cache:
                True
            params:
                length=length,
                temp=temp
            script:
                "src/code/agcrse2_analysis.py"

for seed in range(10):
    for length in [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]:
        rule:
            name:
                f"agcrse2_arrhenius_{length}_{seed}"
            input:
                'src/code/agcrse2_arrhenius.py',
                [f'src/data/reduced/agcrse2/Dmsd{temp}_{length}.h5' for temp in [300, 350, 400, 500, 600, 700]]
            output:
                f'src/data/reduced/agcrse2/bayes_{length}_{seed}.npz'
            conda:
                'environment.yml'
            cache:
                True
            params:
                length=length,
                seed=seed
            script:
                "src/code/agcrse2_arrhenius.py"

rule agcrse2_analysis_generation:
    input:
        'src/scripts/simulation.py',
        'src/scripts/_fig_params.py', 
        [f'src/data/reduced/agcrse2/Dmsd{id}_{length}.h5' for id in [300, 350, 400, 500, 600, 700] for length in [40, 60, 80, 100, 120, 140]],
        [f'src/data/reduced/agcrse2/bayes_{length}_0.npz' for length in [40, 60, 80, 100, 120, 140]]
    output:
        [f'src/tex/figures/agcrse2_{length}.pdf' for length in [40, 60, 80, 100, 120, 140]]
    conda:
        'environment.yml'
    shell:
        "python {input[0]}"

for temp in [500, 600, 700, 800]:
    rule:
        name:
            f"llzo_analysis_{temp}"
        input:
            'src/code/llzo_analysis.py',
            f'src/data/raw/llzo/t{temp}/trajectories{temp}.xyz',
            f'src/data/raw/llzo/t{temp}/traj_{temp}.gro',
            f'src/data/raw/llzo/t{temp}/traj_{temp}.xtc'
        output:
            f'src/data/reduced/llzo/t{temp}/diffusion.h5'
        conda:
            'environment.yml'
        cache:
            True
        params:
            temp=temp
        script:
            "src/code/llzo_analysis.py"

rule llzo_arrhenius:
    input:
        'src/code/llzo_arrhenius.py',
        [f'src/data/reduced/llzo/t{temp}/diffusion.h5' for temp in [500, 600, 700, 800]]
    output:
        'src/data/reduced/llzo/arrhenius_samples.npz'
    conda:
        'environment.yml'
    cache:
        True
    shell:
        "python {input[0]}"

rule make_tex:
    input:
        'src/scripts/make_tex.py',
        ['src/data/reduced/llzo/arrhenius_samples.npz'],
        [f'src/data/reduced/agcrse2/bayes_{i}_{j}.npz' for i in [40, 60, 80, 100, 120, 140] for j in range(10)]
    output:
        'src/tex/output/activation_energy_llzo.txt',
        'src/tex/output/preexp_llzo.txt',
        'src/tex/output/d_300.txt',
        'src/tex/output/bayes_40.txt',
        'src/tex/output/bayes_140.txt',
        [f'src/tex/output/dt{temp}.txt' for temp in [300, 350, 400, 500, 600, 700]]
    conda:
        'environment.yml'
    shell:
        "python {input[0]}"

rule unzip_data:
    input:
        'src/data/raw/agcrse2/agcrse2_T300.zip',
        'src/data/raw/agcrse2/agcrse2_T350.zip',
        'src/data/raw/agcrse2/agcrse2_T400.zip',
        'src/data/raw/agcrse2/agcrse2_T500.zip',
        'src/data/raw/agcrse2/agcrse2_T600.zip',
        'src/data/raw/agcrse2/agcrse2_T700.zip'
    output:
        [f'src/data/raw/agcrse2/agcrse2_T300_{i}.lammpstrj' for i in range(1, 22)],
        [f'src/data/raw/agcrse2/agcrse2_T350_{i}.lammpstrj' for i in range(1, 22)],
        [f'src/data/raw/agcrse2/agcrse2_T400_{i}.lammpstrj' for i in range(1, 22)],
        [f'src/data/raw/agcrse2/agcrse2_T500_{i}.lammpstrj' for i in range(1, 22)],
        [f'src/data/raw/agcrse2/agcrse2_T600_{i}.lammpstrj' for i in range(1, 22)],
        [f'src/data/raw/agcrse2/agcrse2_T700_{i}.lammpstrj' for i in range(1, 22)]
    conda:
        'environment.yml'
    shell:
        "cd src/data/raw/agcrse2 && unzip -o agcrse2_T300.zip && \
         unzip -o agcrse2_T350.zip && \
         unzip -o agcrse2_T400.zip && \
         unzip -o agcrse2_T500.zip && \
         unzip -o agcrse2_T600.zip && \
         unzip -o agcrse2_T700.zip && rm -r __MACOSX && cd ../../../.."
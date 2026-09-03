#methods:
#   - ASHIC: ASHIC-ZIPM model for diploid Hi-C (used in dip3d publication)
#   - LorDG: Lorentzian objective function based reconstruction
#   - HiC-GNN: Graph Convolutional Neural Network based reconstruction
#   - ParticleChromo3D: Particle Swarm Optimization based reconstruction

import os
import glob

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------

def parse_structure3d_roi_config():
    """
    Parse the new structure3d_roi config format and return list of (chr, roi, res) tuples.
    
    Config format:
        structure3d_roi:
          - "chr11"
            - "1000000-3000000" : [50000, 100000, 250000]
          - "chrX":
            - "complete" : [250000]
            - "66919518-77620517" : [50000]
    
    Returns: List of dicts with keys: chr, roi, res
    """
    combinations = []
    roi_config = config.get('structure3d_roi', [])
    
    for chr_entry in roi_config:
        if isinstance(chr_entry, dict):
            for chrom, roi_list in chr_entry.items():
                if isinstance(roi_list, list):
                    for roi_entry in roi_list:
                        if isinstance(roi_entry, dict):
                            for roi, resolutions in roi_entry.items():
                                if isinstance(resolutions, list):
                                    for res in resolutions:
                                        combinations.append({
                                            'chr': chrom,
                                            'roi': roi,
                                            'res': res
                                        })
                                else:
                                    # Single resolution
                                    combinations.append({
                                        'chr': chrom,
                                        'roi': roi,
                                        'res': resolutions
                                    })
                        elif isinstance(roi_entry, str):
                            # ROI without specific resolution, use default
                            default_res = config.get('structure3d_default_resolution', 100000)
                            combinations.append({
                                'chr': chrom,
                                'roi': roi_entry,
                                'res': default_res
                            })
    
    return combinations


# -------------------------------------------------------------------
# Main targets
# -------------------------------------------------------------------

asm_struct3d = ["T2T12", "T2T13"]
#'hicgnn',
tools_struct3d = ['ashic', 'lordg', "particlechromo3d"]


rule all_structure3d:
    input:
        # Run all tools for specified samples, chromosomes, ROIs, and resolutions
        # Combo format: "chr12/1000000-3000000/50000"
        expand("analysis_other/3dstructure/{asm}/{combo}/{tool}/{hp}/structure.pdb",
            asm = asm_struct3d,
            combo = [f"{d['chr']}/{d['roi']}/{d['res']}" for d in parse_structure3d_roi_config()],
            tool = tools_struct3d,
            hp = ['hp1', 'hp2']),

        # 3D Visualizations
        expand("analysis_other/3dstructure/{asm}/{combo}/{tool}/{hp}/visualization.png",
            asm = asm_struct3d,
            combo = [f"{d['chr']}/{d['roi']}/{d['res']}" for d in parse_structure3d_roi_config()],
            tool = tools_struct3d,
            hp = ['hp1', 'hp2']),

        # PyMol Visualizations
        expand("analysis_other/3dstructure/{asm}/{combo}/{tool}/{hp}/pymol_render.png",
            asm = asm_struct3d,
            combo= [f"{d['chr']}/{d['roi']}/{d['res']}" for d in parse_structure3d_roi_config()],
            tool = tools_struct3d,
            hp = ['hp1', 'hp2']),

        # Combined visualizations for every locus
        expand("analysis_other/3dstructure/{asm}/{combo}/combined_visualization.html",
            asm = asm_struct3d,
            combo = [f"{d['chr']}/{d['roi']}/{d['res']}" for d in parse_structure3d_roi_config()]),

        #Combined structure metrics
        "analysis_other/3dstructure/structure3d_summary.tsv"

        



# -------------------------------------------------------------------
# Step 0: Prepare input contact matrices from dip3d haplotagged data
# -------------------------------------------------------------------


# Convert cooler to contact matrix in dense format 
# Does not neeed to be balanced. Uses bin IDs as first two columns (starting at chr1)
rule structure3d_cooler_to_list:
    input:
        cool = "analysis_other/porec/{asm}/cooler/{asm}_{res}_balanced.cool"
    output:
        list="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/contacts/{hp}.contact_list.txt",
    params:
        roi = lambda wc: "" if wc.roi == "complete" else f":{wc.roi}",
    conda:
        "../env/structure3d.yml"
    threads: 4
    log:
        "logs/structure3d/cooler_to_matrix/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/cooler_to_matrix/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        cooler dump \
            -r {wildcards.chr}{params.roi} \
            {input.cool} \
            > {output.list} 2>>{log} \
        """

rule structure3d_list_to_matrix:
    input:
        list = rules.structure3d_cooler_to_list.output.list,
        chrom_sizes = config['ref'] + ".chrom-size.txt"
    output:
        matrix="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/contacts/{hp}.matrix.txt",
    params:
        start=lambda wc: 0 if wc.roi == "complete" else int(wc.roi.split("-")[0]),
        end=lambda wc, input: get_chrom_size(input.chrom_sizes, wc.chr) if wc.roi == "complete" else int(wc.roi.split("-")[1]),
        script = workflow.source_path("../scripts/54_cooler_to_matrix.py")
    conda:
        "../env/structure3d.yml"
    threads: 4
    log:
        "logs/structure3d/cooler_to_matrix/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/cooler_to_matrix/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        cat {input.list} \
        | python {params.script} \
            --res {wildcards.res} \
            --start {params.start} \
            --end {params.end} \
            > {output.matrix} 2>>{log}
        """

# -------------------------------------------------------------------
# ASHIC: ASHIC-ZIPM 3D Structure Prediction
# -------------------------------------------------------------------


rule structure3d_get_chr_sizes:
    input:
        ref=config['ref']
    output:
        txt=config['ref'] + ".chrom-size.txt"
    params:
        dip3d=config['dip3d']
    log:
        "logs/structure3d/ashic_get_chr_sizes.log"
    benchmark:
        "runtimes/structure3d/ashic_get_chr_sizes.txt"
    shell:
        """
        {params.dip3d} extract-ashic-chr-size \
            {input.ref} \
            {output.txt} \
            >{log} 2>&1
        """


rule structure3d_ashic_frag_to_read_pair:
    """Convert haplotagged fragments into ASHIC read-pair format"""
    input:
        frag_list="analysis_other/dip3d/{asm}/4-haplotag/{chr}/imputed-frag-hap-list"
    output:
        pairs="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/ashic_read_pair"
    params:
        dip3d=config['dip3d'],
        min_mapQ=config.get('structure3d_ashic_min_mapq', 5),
        chr=lambda wc: wc.chr
    log:
        "logs/structure3d/ashic_frag_to_read_pair/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/ashic_frag_to_read_pair/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        mkdir -p $(dirname {output.pairs})
        {params.dip3d} frag-to-ashic-read-pair \
            {input.frag_list} \
            {params.chr} \
            {params.min_mapQ} \
            {output.pairs} \
            >{log} 2>&1
        """


rule structure3d_ashic_split2chrs:
    """Split ASHIC read pairs into chromosome-specific bins"""
    input:
        read_pair="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/ashic_read_pair"
    output:
        split_dir=directory("analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/split")
    params:
        ashic_data=config.get('ashic_data', 'bin/ASHIC/ashic/cli/ashic_data.py')
    conda:
        "../env/ashic.yml"
    log:
        "logs/structure3d/ashic_split2chrs/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/ashic_split2chrs/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        mkdir -p {output.split_dir}
        python {params.ashic_data} split2chrs \
            --chr1 1 --allele1 3 \
            --chr2 4 --allele2 6 \
            {input.read_pair} \
            {output.split_dir} \
            >{log} 2>&1
        """


rule structure3d_ashic_binning:
    """Bin ASHIC data for the requested region and resolution"""
    input:
        split_dir="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/split",
        chr_sizes=rules.structure3d_get_chr_sizes.output.txt
    output:
        binned_dir=directory("analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/binned")
    params:
        ashic_data=config.get('ashic_data', 'bin/ASHIC/ashic/cli/ashic_data.py'),
        chr=lambda wc: wc.chr,
        start=lambda wc: "0" if wc.roi == "complete" else wc.roi.split("-")[0],
        end=lambda wc, input: str(get_chrom_size(input.chr_sizes, wc.chr)) if wc.roi == "complete" else wc.roi.split("-")[1],
        res=lambda wc: wc.res
    conda:
        "../env/ashic.yml"
    log:
        "logs/structure3d/ashic_binning/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/ashic_binning/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        python {params.ashic_data} binning \
            --c1=1 --p1=2 --a1=3 \
            --c2=4 --p2=5 --a2=6 \
            --res={params.res} \
            --chrom={params.chr} \
            --start={params.start} \
            --end={params.end} \
            --genome {input.chr_sizes} \
            {input.split_dir}/ashic_read_pair \
            {output.binned_dir} \
            >{log} 2>&1
        """


checkpoint structure3d_ashic_pack:
    """Pack binned ASHIC matrices into pickle format"""
    input:
        binned_dir="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/binned"
    output:
        directory("analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/packed")
    params:
        ashic_data=config.get('ashic_data', 'bin/ASHIC/ashic/cli/ashic_data.py')
    conda:
        "../env/ashic.yml"
    log:
        "logs/structure3d/ashic_pack/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/ashic_pack/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        python {params.ashic_data} pack \
            {input.binned_dir} \
            {output} \
            >{log} 2>&1 || true
        """


def get_structure3d_ashic_pickle_file(wildcards):
    """Resolve the actual ASHIC pickle from the packing checkpoint"""
    import os
    import glob

    checkpoint_output = checkpoints.structure3d_ashic_pack.get(**wildcards).output[0]
    pickle_files = glob.glob(os.path.join(checkpoint_output, "*.pickle"))

    if not pickle_files:
        raise ValueError(f"No pickle file found in {checkpoint_output}")
    if len(pickle_files) > 1:
        pattern = f"ashic_read_pair_{wildcards.chr}_"
        matching = [f for f in pickle_files if pattern in os.path.basename(f)]
        if matching:
            return matching[0]
    return pickle_files[0]


rule structure3d_ashic_run:
    """Run ASHIC-ZIPM to generate diploid 3D structures"""
    input:
        pickle=get_structure3d_ashic_pickle_file
    output:
        html="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/structure_3d.html",
        coords="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/matrices/structure.txt",
    params:
        ashic=config.get('ashic', 'bin/ASHIC/ashic/__main__.py'),
        model=config.get('ashic_model', 'ASHIC-ZIPM')
    conda:
        "../env/ashic.yml"
    threads: 4
    log:
        "logs/structure3d/ashic_run/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/ashic_run/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        python {params.ashic} \
            -i {input.pickle} \
            -o $(dirname {output.html}) \
            --model {params.model} \
            >{log} 2>&1
        """


rule structure3d_ashic_to_pdb:
    """Convert ASHIC output to PDB format"""
    input:
        coords="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/matrices/structure.txt",
        chr_sizes=rules.structure3d_get_chr_sizes.output.txt
    output:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/ashic/{hp}/structure.pdb"
    params:
        chr=lambda wc: wc.chr,
        start=lambda wc: 0 if wc.roi == "complete" else int(wc.roi.split("-")[0]),
        res=lambda wc: int(wc.res)
    log:
        "logs/structure3d/ashic_to_pdb/{asm}.{chr}.{roi}.{res}.{hp}.log"
    run:
        import numpy as np
        
        # Resolve start for complete ROI
        if wildcards.roi == "complete":
            start_pos = 0
        else:
            start_pos = int(wildcards.roi.split("-")[0])
        
        with open(log[0], 'w') as logf:
            try:
                coords = np.loadtxt(input.coords)
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 3)
                
                logf.write(f"Loaded {len(coords)} coordinates\n")
                
                with open(output.pdb, 'w') as pdb:
                    pdb.write(f"HEADER    3D GENOME STRUCTURE - ASHIC\n")
                    pdb.write(f"TITLE     {params.chr}:{start_pos}-{start_pos + len(coords) * params.res}\n")
                    pdb.write(f"REMARK    Resolution: {params.res} bp\n")
                    
                    for i, (x, y, z) in enumerate(coords):
                        # Scale coordinates
                        x, y, z = x * 100, y * 100, z * 100
                        pos = start_pos + i * params.res
                        pdb.write(f"ATOM  {i+1:5d}  CA  ALA A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                    
                    # Add CONECT records for chain connectivity
                    for i in range(len(coords) - 1):
                        pdb.write(f"CONECT{i+1:5d}{i+2:5d}\n")
                    
                    pdb.write("END\n")
                    
            except Exception as e:
                logf.write(f"Error: {e}\n")
                # Create minimal valid PDB
                with open(output.pdb, 'w') as pdb:
                    pdb.write("HEADER    3D GENOME STRUCTURE - ASHIC (EMPTY)\n")
                    pdb.write("END\n")


# -------------------------------------------------------------------
# LorDG: Lorentzian Objective Function 3D Reconstruction
# -------------------------------------------------------------------

rule structure3d_lordg_prepare:
    """Prepare LorDG parameter file"""
    input:
        contact_list=rules.structure3d_cooler_to_list.output.list
    output:
        params="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/lordg/{hp}/parameters.txt"
    params:
        output_folder=lambda wc, output: os.path.dirname(output.params),
        num_models=config.get('lordg_num_models', 1),
        learning_rate=config.get('lordg_learning_rate', 1),
        max_iteration=config.get('lordg_max_iteration', 10000)
    log:
        "logs/structure3d/lordg_prepare/{asm}.{chr}.{roi}.{res}.{hp}.log"
    run:
        import os
        
        with open(log[0], 'w') as logf:
            os.makedirs(params.output_folder, exist_ok=True)
            logf.write(f"Creating LorDG parameters file\n")
            
            with open(output.params, 'w') as f:
                f.write(f"NUM = {params.num_models}\n")
                f.write(f"OUTPUT_FOLDER = {params.output_folder}/output\n")
                f.write(f"INPUT_FILE = {input.contact_list}\n")
                f.write(f"VERBOSE = true\n")
                f.write(f"LEARNING_RATE = {params.learning_rate}\n")
                f.write(f"MAX_ITERATION = {params.max_iteration}\n")
            
            logf.write(f"Parameters written to {output.params}\n")


rule structure3d_lordg_run:
    """
    Run LorDG 3D structure prediction
    Input format: contact list with three numeric columns:
    position_1 position_2 interaction_frequencies
    """

    input:
        params="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/lordg/{hp}/parameters.txt",
        contact_list="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/contacts/{hp}.contact_list.txt"
    output:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/lordg/{hp}/structure.pdb",
        log_out="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/lordg/{hp}/lordg_log.txt"
    params:
        lordg_jar=config.get('lordg_jar', 'bin/LorDG/bin/3DDistanceBaseLorentz.jar'),
        output_folder=lambda wc, output: os.path.dirname(output.pdb) + "/output"
    conda:
        "../env/lordg.yml"
    threads: 4
    log:
        "logs/structure3d/lordg_run/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/lordg_run/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        mkdir -p {params.output_folder}
        
        # Run LorDG
        java -jar {params.lordg_jar} {input.params} >{log} 2>&1

        # Copy first generated PDB and log file
        cp $(ls {params.output_folder}/*.pdb | head -1) {output.pdb} 2>>{log}
        cp $(ls {params.output_folder}/*log*.txt 2>/dev/null | head -1) {output.log_out} 2>>{log} || echo "No log file found" > {output.log_out}
       
        """


# -------------------------------------------------------------------
# HiC-GNN: Graph Convolutional Neural Network 3D Reconstruction  
# -------------------------------------------------------------------

rule structure3d_hicgnn:
    input:
        matrix="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/contacts/{hp}.contact_list.txt"
    output:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/hicgnn/{hp}/structure.pdb"
    resources:
        hicgnn_slot=1
    params:
        hicgnn_docker=config.get('hicgnn_docker', 'oluwadarelab/hicgnn:latest'),
        hicgnn_repo=config.get('hicgnn_repo', 'bin/HiC-GNN'),
        conversions=config.get('hicgnn_conversions', '[.1,.1,2]'),
        epochs=config.get('hicgnn_epochs', 10),
        learning_rate=config.get('hicgnn_learning_rate', 0.001)
    conda:
        "../env/hicgnn.yml"
    threads: 8
    log:
        "logs/structure3d/hicgnn/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/hicgnn/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        WD=$(realpath $(dirname {output.pdb}))
        mkdir -p $WD/output

        # Copy input matrix to work directory
        cp {input.matrix} $WD/input_matrix.txt
        
        # Run HiC-GNN using Docker
        docker run --rm \
            --user $(id -u):$(id -g) \
            -v $(realpath {params.hicgnn_repo}):/HiC-GNN \
            -v $WD:/workdir \
            {params.hicgnn_docker} \
            bash -c "
                python HiC-GNN_main.py \
                    /workdir/input_matrix.txt \
                    -c '{params.conversions}' \
                    -ep {params.epochs} \
                    -lr {params.learning_rate} \
                && mv -v Outputs/* /workdir/output/ \
                && rm -v Data/input_matrix* 
            " >{log} 2>&1
        
        # Copy output PDB
        cp $WD/output/input_matrix_structure.pdb {output.pdb} 2>>{log} 
        """

# -------------------------------------------------------------------
# ParticleChromo3D: Particle Swarm Optimization-based 3D Reconstruction
# -------------------------------------------------------------------

rule structure3d_particlechromo3d_run:
    input:
        matrix=rules.structure3d_list_to_matrix.output.matrix
    output:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/particlechromo3d/{hp}/structure.pdb",
    params:
        particlechromo3d=config.get('particlechromo3d_path', 'ParticleChromo3D'),
        swarm_size=config.get('particlechromo3d_swarm_size', 15),
        iterations=config.get('particlechromo3d_iterations', 30000),
        threshold=config.get('particlechromo3d_threshold', 0.000001),
        randrange=config.get('particlechromo3d_randrange', 1),
        loss_function=config.get('particlechromo3d_loss_function', 2),  # RMSE
    conda:
        "../env/particlechromo3d.yml"
    threads: 4
    log:
        "logs/structure3d/particlechromo3d/{asm}.{chr}.{roi}.{res}.{hp}.log"
    benchmark:
        "runtimes/structure3d/particlechromo3d/{asm}.{chr}.{roi}.{res}.{hp}.txt"
    shell:
        """
        # ParticleChromo3D expects normalized contact matrix
        python {params.particlechromo3d} \
            -ss {params.swarm_size} \
            -itt {params.iterations} \
            -t {params.threshold} \
            -rr {params.randrange} \
            -o $(dirname {output.pdb})/$(basename {output.pdb} .pdb) \
            -lf {params.loss_function} \
            {input.matrix} \
            >{log} 2>&1
        """


# -------------------------------------------------------------------
# Visualization and Plotting
# -------------------------------------------------------------------

rule structure3d_visualize:
    """Create 3D visualization for a single structure"""
    input:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/structure.pdb"
    output:
        png="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/visualization.png",
        html="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/visualization.html"
    params:
        title=lambda wc: f"{wc.asm} {wc.chr}:{wc.roi} {wc.res}bp {wc.tool} {wc.hp}"
    conda:
        "../env/structure3d_viz.yml"
    log:
        "logs/structure3d/visualize/{asm}.{chr}.{roi}.{res}.{tool}.{hp}.log"
    shell:
        """
        python3 workflow/scripts/49_visualize_3d_structure.py \
            --pdb {input.pdb} \
            --output-png {output.png} \
            --output-html {output.html} \
            --title "{params.title}" \
            >{log} 2>&1
        """


rule structure3d_combined_visualization:
    """Create combined visualization comparing all tools and haplotypes"""
    input:
        pdbs=lambda wc: expand("analysis_other/3dstructure/{{asm}}/{{chr}}/{{roi}}/{{res}}/{tool}/{hp}/structure.pdb",
            tool=tools_struct3d,
            hp=['hp1', 'hp2'])
    output:
        html="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/combined_visualization.html",
        comparison="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/method_comparison.tsv"
    params:
        title=lambda wc: f"{wc.asm} {wc.chr}:{wc.roi} {wc.res}bp - Method Comparison"
    conda:
        "../env/structure3d_viz.yml"
    log:
        "logs/structure3d/combined_viz/{asm}.{chr}.{roi}.{res}.log"
    shell:
        """
        python3 workflow/scripts/50_compare_3d_structures.py \
            --pdbs {input.pdbs} \
            --output-html {output.html} \
            --output-comparison {output.comparison} \
            --title "{params.title}" \
            >{log} 2>&1
        """


rule structure3d_pymol_render:
    """Optional: High-quality PyMOL rendering"""
    input:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/structure.pdb"
    output:
        png="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/pymol_render.png"
    conda:
        "../env/pymol.yml"
    log:
        "logs/structure3d/pymol_render/{asm}.{chr}.{roi}.{res}.{tool}.{hp}.log"
    shell:
        """
        # Create PyMOL script
        cat > $(dirname {output.png})/render.pml << 'PYMOLSCRIPT'
load {input.pdb}
bg_color white
set max_threads, 1
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1
spectrum count, rainbow_rev, all
set cartoon_tube_radius, 1.0
show cartoon
hide lines
set ray_shadow, 0
ray 2000, 2000
png {output.png}, dpi=300
quit
PYMOLSCRIPT

        # Run PyMOL in headless mode
        pymol -cq $(dirname {output.png})/render.pml >{log} 2>&1 || true
    
        """


# -------------------------------------------------------------------
# Quality metrics and comparison
# -------------------------------------------------------------------

rule structure3d_calculate_metrics:
    """Calculate quality metrics for 3D structure"""
    input:
        pdb="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/structure.pdb",
        matrix="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/contacts/{hp}.matrix.txt"
    output:
        metrics="analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/metrics.json"
    conda:
        "../env/structure3d.yml"
    log:
        "logs/structure3d/metrics/{asm}.{chr}.{roi}.{res}.{tool}.{hp}.log"
    shell:
        """
        python3 workflow/scripts/51_calculate_structure_metrics.py \
            --pdb {input.pdb} \
            --contact-matrix {input.matrix} \
            --output {output.metrics} \
            >{log} 2>&1
        """


rule structure3d_summary:
    """Create summary of all 3D structure predictions"""
    input:
        metrics=expand("analysis_other/3dstructure/{asm}/{combo}/{tool}/{hp}/metrics.json",
            asm = asm_struct3d,
            combo = [f"{d['chr']}/{d['roi']}/{d['res']}" for d in parse_structure3d_roi_config()],
            tool = tools_struct3d,
            hp = ['hp1', 'hp2']),
    output:
        summary="analysis_other/3dstructure/structure3d_summary.tsv"
    conda:
        "../env/structure3d.yml"
    log:
        "logs/structure3d/summary.log"
    shell:
        """
        python3 workflow/scripts/52_summarize_structure3d.py \
            --metrics {input.metrics} \
            --output {output.summary} \
            >{log} 2>&1
        """

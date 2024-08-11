importScripts("https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js");

function sendPatch(patch, buffers, msg_id) {
  self.postMessage({
    type: 'patch',
    patch: patch,
    buffers: buffers
  })
}

async function startApplication() {
  console.log("Loading pyodide!");
  self.postMessage({type: 'status', msg: 'Loading pyodide'})
  self.pyodide = await loadPyodide();
  self.pyodide.globals.set("sendPatch", sendPatch);
  console.log("Loaded!");
  await self.pyodide.loadPackage("micropip");
  const env_spec = ['https://cdn.holoviz.org/panel/wheels/bokeh-3.4.3-py3-none-any.whl', 'https://cdn.holoviz.org/panel/1.4.5/dist/wheels/panel-1.4.5-py3-none-any.whl', 'pyodide-http==0.2.1', 'matplotlib', 'numpy', 'pandas']
  for (const pkg of env_spec) {
    let pkg_name;
    if (pkg.endsWith('.whl')) {
      pkg_name = pkg.split('/').slice(-1)[0].split('-')[0]
    } else {
      pkg_name = pkg
    }
    self.postMessage({type: 'status', msg: `Installing ${pkg_name}`})
    try {
      await self.pyodide.runPythonAsync(`
        import micropip
        await micropip.install('${pkg}');
      `);
    } catch(e) {
      console.log(e)
      self.postMessage({
	type: 'status',
	msg: `Error while installing ${pkg_name}`
      });
    }
  }
  console.log("Packages loaded!");
  self.postMessage({type: 'status', msg: 'Executing code'})
  const code = `
  \nimport asyncio\n\nfrom panel.io.pyodide import init_doc, write_doc\n\ninit_doc()\n\nfrom pathlib import Path\nfrom collections import defaultdict\nimport numpy as np\nfrom time import time\nfrom multiprocessing import Process, Array\nfrom math import ceil\nimport matplotlib.pyplot as plt\nimport pandas as pd\nimport panel as pn\nfrom io import StringIO\nfrom datetime import datetime\n\n\n\npn.extension()\npn.extension('tabulator')\n\n\n\n# Globals\n\ngeno_file_path = ''\n\npopulations_dict = defaultdict(list)\npopulations = []\n\nparsed_sel_pops = []\nsel_pops = []\nsnp_names = []\n\nnum_alleles = 0\nallele_frequencies = defaultdict(list)\n\ninvalid_indices = []\nnum_valid_alleles = 0\n\nroles_pops = []\n\nhybrid = ''\nparent1 = ''\nparent2 = ''\nauxiliaries = []\n\nalpha_pre_jl = None\n\nf4ab_prime = None\nf4xb_prime = None\nalpha = None\nalpha_error = None\n\nf4ab_std = None\nf4xb_std = None\nalpha_std = None\nalpha_std_error = None\n\nalpha_ratio = None\nalpha_ratio_avg = None\nalpha_ratio_std_dev = None\nalpha_ratio_hist = None\n\ncosine_pre_jl = None\nangle_pre_jl = None\npercentage_pre_jl = None\n\ncosine_post_jl = None\nangle_post_jl = None\npercentage_post_jl = None\n\nnum_cases = None\n\nf3_test = None\n\n\n\n# Input and output widgets\n\nfile_selector = pn.widgets.FileSelector(name = 'Select .geno, .ind and .snp files', only_files = True)\nload_files_button = pn.widgets.Button(name = 'Parse and check input files', button_type = 'primary')\n\nalert_pane = pn.pane.Alert('### Input files selection\\nPlease, select a triad of .geno, .ind and .snp input files, and optionally a .dat input file with a list of selected populations.\\nThen press the parse button.', alert_type = 'primary')\n\nreset_sel_pops_button = pn.widgets.Button(name = 'Reset selected populations', button_type = 'primary')\n\nnum_procs_input = pn.widgets.IntInput(name = 'Number of parallel processes', value = 1, step = 1, start = 1)\ncompute_freqs_button = pn.widgets.Button(name = 'Compute allele frequencies', button_type = 'primary', disabled = True)\nfreqs_download_button = pn.widgets.FileDownload(label = 'Download allele frequencies', disabled = True, button_type = 'primary')\n\nhybrid_select = pn.widgets.Select(name = 'Hybrid', options = [], size = 10)\nparent1_select = pn.widgets.Select(name = 'Parent 1', options = [], size = 10)\nparent2_select = pn.widgets.Select(name = 'Parent 2', options = [], size = 10)\naux_pops_select = pn.widgets.MultiSelect(name = 'Auxiliaries', options = [], size = 10)\n\n\nplot_width_input = pn.widgets.FloatInput(name = 'Plot width (inches)', value = 4, step = 0.1, start = 0.1)\nplot_height_input = pn.widgets.FloatInput(name = 'Plot height (inches)', value = 4, step = 0.1, start = 0.1)\nplot_title_size_input = pn.widgets.IntInput(name = 'Plot title font size', value = 10, step = 1, start = 1)\nplot_label_size_input = pn.widgets.IntInput(name = 'Plot labels font size', value = 10, step = 1, start = 1)\ncompute_results_button = pn.widgets.Button(name = 'Compute results', button_type = 'primary', disabled = True)\nresults_download_button = pn.widgets.FileDownload(label = 'Download results', disabled = True, button_type = 'primary')\nf4_points_download_button = pn.widgets.FileDownload(label = 'Download f4-points', disabled = True, button_type = 'primary')\nresulting_data_output = pn.pane.Markdown()\n\n\n\n# Tables\n\navail_pops_filter = pn.widgets.TextInput(name = 'Search populations:', placeholder = 'Enter population name', disabled = True)\navail_pops_table = pn.widgets.Tabulator(show_index = False, disabled = True, selectable = False, pagination = 'local', page_size = 10, align = ('start'), widths = {'population': '100%'})\nsel_pops_table = pn.widgets.Tabulator(show_index = False, disabled = True, selectable = False, sortable = False, pagination = 'local', page_size = 10, align = ('start'), widths = {'population': '100%'}, buttons = {'remove': "<i class='fa fa-times'></i>"})\n\n\n\n# Plot panes\n\npane_margin = 10\n\nf4prime_fit_pane = pn.pane.Matplotlib(align = ('center', 'start'), margin = pane_margin)\nf4_fit_pane = pn.pane.Matplotlib(align = ('center', 'start'), margin = pane_margin)\nf4_ratio_histogram_pane = pn.pane.Matplotlib(align = ('center', 'start'), margin = pane_margin)\n\n\n\ndef set_alert_pane(text, type):\n    alert_pane.object = text\n    alert_pane.alert_type = type\n\n\n\ndef reset_alert_pane(event):\n    text = '### Input files selection\\nPlease, select a triad of .geno, .ind and .snp input files, and optionally a .dat input file with a list of selected populations.\\nThen press the parse button.'\n    alert_pane.object = text\n    alert_pane.alert_type = 'primary'\n\n\n\nfile_selector.param.watch(reset_alert_pane, 'value')\n\n\n\ndef set_pops_table_data(pops, table):\n    idx = list(range(len(pops)))\n    df = pd.DataFrame({'population': pops}, index = idx)\n    table.value = df\n\n\n\ndef contains_filter(df, pattern, column):\n    if not pattern:\n        return df\n    return df[df[column].str.contains(pattern, case = False)]\n\n\n\navail_pops_table.add_filter(pn.bind(contains_filter, pattern = avail_pops_filter, column = 'population'))\n\n\n\ndef set_sel_pops_data(event):\n    global sel_pops\n\n    pop = event.value\n    if pop in sel_pops:\n        sel_pops.remove(pop)\n    else:\n        sel_pops.insert(0, pop)\n\n    set_pops_table_data(sel_pops, sel_pops_table)\n\n    if len(sel_pops) > 0:\n        compute_freqs_button.disabled = False\n    else:\n        compute_freqs_button.disabled = True\n\n\n\navail_pops_table.on_click(set_sel_pops_data)\n\n\n\ndef remove_sel_pop(event):\n    if event.column != 'remove':\n        return\n\n    global sel_pops\n    index = int(event.row)\n    sel_pops.pop(index)\n\n    set_pops_table_data(sel_pops, sel_pops_table)\n\n    if len(sel_pops) == 0:\n        compute_freqs_button.disabled = True\n\n\n\nsel_pops_table.on_click(remove_sel_pop)\n\n\n\ndef set_admixture_model():\n    global hybrid, parent1, parent2, auxiliaries\n\n    hybrid = hybrid_select.value\n    parent1 = parent1_select.value\n    parent2 = parent2_select.value\n    auxiliaries = aux_pops_select.value\n\n\n\ndef init_selects_options():\n    global roles_pops\n    roles_pops = [pop for pop in sel_pops]\n\n    hybrid_select.options = roles_pops\n    parent1_select.options = roles_pops\n    parent2_select.options = roles_pops\n    aux_pops_select.options = roles_pops\n\n    hybrid_select.value = roles_pops[0]\n    parent1_select.value = roles_pops[1]\n    parent2_select.value = roles_pops[2]\n    aux_pops_select.value = roles_pops[3:]\n\n    set_admixture_model()\n\n\n\ndef set_selects_values(event):\n    select_name = event.obj.name\n    old_pop = event.old\n    new_pop = event.new\n\n    if select_name == 'Auxiliaries':\n        pops = [hybrid_select.value, parent1_select.value, parent2_select.value]\n        aux_pops_select.value = [pop for pop in new_pop if pop not in pops]\n    else:\n        if new_pop in aux_pops_select.value:\n            aux_pops_select.value = [pop for pop in aux_pops_select.value if pop != new_pop]\n        else:\n            selects = []\n            if select_name == 'Hybrid':\n                selects = [parent1_select, parent2_select]\n            elif select_name == 'Parent 1':\n                selects = [hybrid_select, parent2_select]\n            elif select_name == 'Parent 2':\n                selects = [hybrid_select, parent1_select]\n\n            for sel in selects:\n                if sel.value == new_pop and old_pop != None:\n                    sel.value = old_pop\n\n    set_admixture_model()\n\n    compute_results_button.disabled = (len(aux_pops_select.value) == 0)\n\n\n\nhybrid_select.param.watch(set_selects_values, 'value')\nparent1_select.param.watch(set_selects_values, 'value')\nparent2_select.param.watch(set_selects_values, 'value')\naux_pops_select.param.watch(set_selects_values, 'value')\n\n\n\ndef parse_selected_populations(file_path):\n    global parsed_sel_pops\n    read_sel_pops = []\n\n    with file_path.open(mode = 'r', encoding = 'utf-8') as file:\n        content = file.read()\n        pop_lines = [line.split() for line in content.splitlines() if not line.startswith('#')]\n        read_sel_pops = [pop_line[0] for pop_line in pop_lines]\n\n    invalid_pops = []\n    for pop in read_sel_pops:\n        if not pop in populations:\n            invalid_pops.append(pop)\n\n    global sel_pops\n    sel_pops = [pop for pop in read_sel_pops if pop not in invalid_pops]\n\n    global parsed_sel_pops\n    parsed_sel_pops = [pop for pop in read_sel_pops if pop not in invalid_pops]\n\n    if len(sel_pops) > 0:\n        set_pops_table_data(sel_pops, sel_pops_table)\n        compute_freqs_button.disabled = False\n    else:\n        compute_freqs_button.disabled = True\n\n\n\ndef reset_sel_pops(event):\n    global sel_pops\n    sel_pops = [pop for pop in parsed_sel_pops]\n    set_pops_table_data(sel_pops, sel_pops_table)\n\n    if len(sel_pops) == 0:\n        compute_freqs_button.disabled = True\n        freqs_download_button.disabled = True\n\n\n\npn.bind(reset_sel_pops, reset_sel_pops_button, watch = True)\n\n\n\ndef parse_populations(file_path):\n    global populations_dict, populations\n\n    populations_dict = defaultdict(list)\n    populations = []\n\n    num_rows = 0\n\n    text = f'### Parsing and checking {file_path}\\n'\n    text_lines = [text, '']\n\n    with file_path.open(mode = 'r', encoding = 'utf-8') as file:\n        content = file.read()\n        for index, line in enumerate(content.splitlines()):\n            columns = line.split()\n            pop_name = columns[-1]\n            populations_dict[pop_name].append(index)\n\n            num_rows += 1\n            if (num_rows % 1000 == 0):\n                text_lines[-1] = f'Number of rows: {num_rows}'\n                set_alert_pane('\\n'.join(text_lines), 'warning')\n\n        populations = list(populations_dict.keys())\n\n    return num_rows\n\n\n\ndef parse_snp_names(file_path):\n    global snp_names\n    snp_names = []\n\n    num_alleles = 0\n\n    text = f'### Parsing and checking {file_path}\\n'\n    text_lines = [text, '']\n\n    with file_path.open(mode = 'r', encoding = 'utf-8') as file:\n        content = file.read()\n        for line in content.splitlines():\n            columns = line.split()\n            snp_names.append(columns[0])\n\n            num_alleles += 1\n            if (num_alleles % 10000 == 0):\n                text_lines[-1] = f'Number of alleles: {num_alleles}'\n\n                pn.bind(load_input_files, load_files_button, watch = True)\n\n                set_alert_pane('\\n'.join(text_lines), 'warning')\n\n    return num_alleles\n\n\n\ndef geno_table_shape(file_path):\n    num_rows = 0\n    num_columns = []\n\n    text = f'### Parsing and checking {file_path}\\n'\n    text_lines = [text, '']\n\n    with file_path.open(mode = 'r', encoding = 'utf-8') as file:\n        for row in file:\n            num_rows += 1\n            num_columns.append(len(row) - 1)\n\n            if (num_rows % 1000 == 0):\n                text_lines[-1] = f'Number of rows: {num_rows}'\n                set_alert_pane('\\n'.join(text_lines), 'warning')\n\n        text_lines[-1] = f'Number of rows: {num_rows}'\n        set_alert_pane('\\n'.join(text_lines), 'warning')\n\n    return num_rows, num_columns\n\n\n\ndef load_input_files(event):\n    file_paths = file_selector.value\n\n    # Check selected files for invalid suffixes\n\n    invalid_file_paths = []\n\n    for fp in file_paths:\n        file_path = Path(fp)\n        if file_path.suffix not in ['.geno', '.ind', '.snp', '.dat']:\n            invalid_file_paths.append(fp + '\\n')\n\n    if len(invalid_file_paths) > 0:\n        invalid_fp_list = '- ' + '- '.join(invalid_file_paths)\n        text = f'### Unrecognized file suffixes\\nValid input files end with the suffixes: .geno, .ind and .snp\\nPlease, deselect the following files:\\n{invalid_fp_list}'\n        set_alert_pane(text, 'danger')\n        return\n\n    # Check for triad\n\n    suffixes = [Path(fp).suffix for fp in file_paths]\n\n    count_geno = suffixes.count('.geno')\n    count_ind = suffixes.count('.ind')\n    count_snp = suffixes.count('.snp')\n    count_dat = suffixes.count('.dat')\n\n    if count_geno != 1 or count_ind != 1 or count_snp != 1 or count_dat > 1:\n        text = '### Wrong number of input files\\nPlease, select a triad of input files (.geno, .ind, .snp) and optionally a single selected populations input file (.dat)'\n        set_alert_pane(text, 'danger')\n        return\n\n    # Succesful selection\n\n    global geno_file_path\n\n    geno_file_path = next((Path(fp) for fp in file_paths if Path(fp).suffix == '.geno'), None)\n    ind_file_path = next((Path(fp) for fp in file_paths if Path(fp).suffix == '.ind'), None)\n    snp_file_path = next((Path(fp) for fp in file_paths if Path(fp).suffix == '.snp'), None)\n    dat_file_path = next((Path(fp) for fp in file_paths if Path(fp).suffix == '.dat'), None)\n\n    # Parsing\n\n    # Check .geno table\n\n    num_geno_rows, num_geno_columns = geno_table_shape(geno_file_path)\n\n    if not all(nc == num_geno_columns[0] for nc in num_geno_columns):\n        text = f'### Parsing failed\\nIn {geno_file_path}: Not all rows are of equal number of columns'\n        set_alert_pane(text, 'danger')\n        return\n\n    # Parse and check number of alleles and rows\n\n    global num_alleles\n\n    num_alleles = parse_snp_names(snp_file_path)\n\n    if num_alleles != num_geno_rows:\n        text = f'### Parsing failed\\nNumber of alleles ({num_alleles}) in .snp file is not equal to number of rows ({num_geno_rows}) in .geno file'\n        set_alert_pane(text, 'danger')\n        return\n\n    # Parse and check columns\n\n    num_ind_rows = parse_populations(ind_file_path)\n\n    if num_ind_rows != num_geno_columns[0]:\n        text = f'### Parsing failed\\nNumber of rows ({num_ind_rows}) in .ind file is not equal to number of columns ({num_geno_columns[0]}) in .geno file'\n        set_alert_pane(text, 'danger')\n        return\n\n    # Parse and check selected populations\n\n    global sel_pops\n    global parsed_sel_pops\n\n    if dat_file_path is not None:\n        parse_selected_populations(dat_file_path)\n    else:\n        sel_pops = []\n        parsed_sel_pops = []\n        compute_freqs_button.disabled = True\n\n    freqs_download_button.disabled = True\n\n    # Set tables\n\n    set_pops_table_data(populations, avail_pops_table)\n    set_pops_table_data(sel_pops, sel_pops_table)\n\n    avail_pops_filter.disabled = False\n\n    text = f'### Parsing successful\\nParsed input files seem to have a correct structure:\\n- {num_geno_rows} rows and {num_geno_columns[0]} columns in {geno_file_path}\\n- {num_alleles} alleles in {snp_file_path}\\n- {num_ind_rows} rows and {len(populations_dict)} populations in {ind_file_path}'\n    if len(sel_pops) > 0:\n        text += f'\\n- {len(sel_pops)} selected populations in {dat_file_path}'\n    set_alert_pane(text, 'success')\n\n\n\npn.bind(load_input_files, load_files_button, watch = True)\n\n\n\ndef invalid_allele_indices():\n    indices = []\n\n    for pop, freqs in allele_frequencies.items():\n        for index, freq in enumerate(freqs):\n            if freq == -1:\n                indices.append(index)\n\n    global invalid_indices\n    invalid_indices = np.unique(np.array(indices, dtype = int))\n\n    global num_valid_alleles\n    num_valid_alleles = num_alleles - invalid_indices.size\n\n\n\ndef remove_invalid_alleles():\n    global allele_frequencies\n    for pop, freqs in allele_frequencies.items():\n        allele_frequencies[pop] = np.delete(freqs, invalid_indices)\n\n\n\ndef allele_frequency(alleles):\n    freq = 0\n    num_alleles = 0\n\n    for a in alleles:\n        if a != 9:\n            freq += (2 - a) / 2\n            num_alleles += 1\n\n    if num_alleles == 0:\n        return -1\n\n    return freq / num_alleles\n\n\n\ndef population_allele_frequencies(pop_indices, allele_freqs):\n    with geno_file_path.open(mode = 'r', encoding = 'utf-8') as file:\n        for index, row in enumerate(file):\n            allele_freqs[index] = allele_frequency([int(row[i]) for i in pop_indices])\n\n\n\ndef parallel_compute_populations_frequencies(event):\n    pop_indices = [populations_dict[pop] for pop in sel_pops]\n    num_computations = len(pop_indices)\n    num_procs = num_procs_input.value\n    batch_size = ceil(num_computations / num_procs)\n    index = 0\n\n    text_lines = ['### Computing']\n    text = f'Computing {num_alleles} frequencies per population for {num_computations} populations in {batch_size} batches of {num_procs} parallel processes.'\n    text_lines.append(text)\n    set_alert_pane('\\n'.join(text_lines), 'warning')\n\n    pops_names = ''\n    text_lines.append('')\n\n    allele_freqs = [Array('d', num_alleles) for i in range(num_computations)]\n\n    t1 = time()\n\n    for nb in range(batch_size):\n        procs = []\n        computing_pops = []\n\n        for nt in range(num_procs):\n            if index < num_computations:\n                p = Process(target = population_allele_frequencies, args = (pop_indices[index], allele_freqs[index]))\n                procs.append(p)\n                p.start()\n                computing_pops.append(sel_pops[index])\n            index += 1\n\n        pops_names += '{' + ', '.join(computing_pops) + '} '\n        text_lines[-1] = pops_names\n        set_alert_pane('\\n'.join(text_lines), 'warning')\n\n        for p in procs:\n            p.join()\n\n    global allele_frequencies\n    allele_frequencies = defaultdict(list)\n    for i, freqs in enumerate(allele_freqs):\n        allele_frequencies[sel_pops[i]] = np.array(freqs)\n\n    t2 = time()\n    comp_time = t2 - t1\n\n    text_lines[0] = '### Finished'\n    text = f'Computation took {comp_time:.2f} seconds.'\n    text_lines.append(text)\n    set_alert_pane('\\n'.join(text_lines), 'success')\n\n    invalid_allele_indices()\n\n    if len(invalid_indices) > 0:\n        text = f'Number of excluded SNPs: {len(invalid_indices)}'\n        text_lines.append(text)\n        set_alert_pane('\\n'.join(text_lines), 'success')\n\n        remove_invalid_alleles()\n\n    init_selects_options()\n\n    compute_results_button.disabled = False\n    resulting_data_output.object = ''\n\n    f4prime_fit_pane.object = None\n    f4_fit_pane.object = None\n    f4_ratio_histogram_pane.object = None\n\n    freqs_download_button.disabled = False\n    results_download_button.disabled = True\n    f4_points_download_button.disabled = True\n\n\n\n\npn.bind(parallel_compute_populations_frequencies, compute_freqs_button, watch = True)\n\n\n\ndef mixing_coefficient_pre_jl(hybrid_freqs, parent1_freqs, parent2_freqs):\n    parent_diff = parent1_freqs - parent2_freqs\n    return np.dot(hybrid_freqs - parent2_freqs, parent_diff) / np.dot(parent_diff, parent_diff)\n\n\n\ndef admixture_angle_pre_jl(hybrid_freqs, parent1_freqs, parent2_freqs):\n    xa = hybrid_freqs - parent1_freqs\n    xb = hybrid_freqs - parent2_freqs\n\n    cosine = np.dot(xa, xb) / np.sqrt(np.dot(xa, xa) * np.dot(xb, xb))\n    angle = np.arccos(cosine)\n    percentage = angle / np.pi\n\n    return cosine, angle * 180 / np.pi, percentage\n\n\n\ndef f3(hybrid_freqs, parent1_freqs, parent2_freqs):\n    num_alleles = hybrid_freqs.size\n    return np.dot(hybrid_freqs - parent1_freqs, hybrid_freqs - parent2_freqs) / num_alleles\n\n\n\ndef f4_prime(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):\n    num_aux_pops = len(aux_freqs)\n    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)\n\n    f4ab_prime = np.zeros(num_pairs)\n    f4xb_prime = np.zeros(num_pairs)\n\n    ab = parent1_freqs - parent2_freqs\n    xb = hybrid_freqs - parent2_freqs\n\n    index = 0\n\n    for i in range(num_aux_pops):\n        for j in range(i + 1, num_aux_pops):\n            ij = aux_freqs[i] - aux_freqs[j]\n            norm_ij = np.linalg.norm(ij)\n            f4ab_prime[index] = np.dot(ab, ij) / norm_ij\n            f4xb_prime[index] = np.dot(xb, ij) / norm_ij\n            index += 1\n\n    return f4ab_prime, f4xb_prime\n\n\n\ndef f4_std(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):\n    num_aux_pops = len(aux_freqs)\n    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)\n\n    f4ab_std = np.zeros(num_pairs)\n    f4xb_std = np.zeros(num_pairs)\n\n    ab = parent1_freqs - parent2_freqs\n    xb = hybrid_freqs - parent2_freqs\n\n    index = 0\n\n    for i in range(num_aux_pops):\n        for j in range(i + 1, num_aux_pops):\n            ij = aux_freqs[i] - aux_freqs[j]\n            f4ab_std[index] = np.dot(ab, ij)\n            f4xb_std[index] = np.dot(xb, ij)\n            index += 1\n\n    num_alleles = hybrid_freqs.size\n\n    return f4ab_std / num_alleles, f4xb_std / num_alleles\n\n\n\ndef least_squares(x, y):\n    dim = len(x)\n\n    A = np.vstack([x, np.zeros(dim)]).T\n    alpha = np.linalg.lstsq(A, y)[0][0]\n\n    Q = 0\n    for i in range(dim):\n        Q += (y[i] - alpha * x[i]) ** 2\n\n    x_avg = 0\n    for i in range(dim):\n        x_avg += x[i]\n    x_avg /= dim\n\n    x_dev = 0\n    for i in range(dim):\n        x_dev += (x[i] - x_avg) ** 2\n\n    s_alpha = np.sqrt(Q / ((dim - 2) * x_dev))\n    t = 1.98\n\n    error = s_alpha * t\n\n    return alpha, error\n\n\n\ndef admixture_angle_post_jl(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):\n    num_aux_pops = len(aux_freqs)\n    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)\n\n    xa = hybrid_freqs - parent1_freqs\n    xb = hybrid_freqs - parent2_freqs\n\n    sum1 = 0\n    sum2 = 0\n    sum3 = 0\n\n\n    for i in range(num_aux_pops):\n        for j in range(i + 1, num_aux_pops):\n            ij = aux_freqs[i] - aux_freqs[j]\n\n            xaij = np.dot(xa, ij)\n            xbij = np.dot(xb, ij)\n            ijij = np.dot(ij, ij)\n\n            sum1 += xaij * xbij / ijij\n            sum2 += (xaij ** 2) / ijij\n            sum3 += (xbij ** 2) / ijij\n\n    cosine = sum1 / np.sqrt(sum2 * sum3)\n    angle = np.arccos(cosine)\n    percentage = angle / np.pi\n\n    return cosine, angle * 180 / np.pi, percentage\n\n\n\ndef f4_ratio(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):\n    num_aux_pops = len(aux_freqs)\n    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)\n\n    xb = hybrid_freqs - parent2_freqs\n    ab = parent1_freqs - parent2_freqs\n\n    alpha = np.zeros(num_pairs)\n\n    index = 0\n\n    for i in range(num_aux_pops):\n        for j in range(i + 1, num_aux_pops):\n            ij = aux_freqs[i] - aux_freqs[j]\n            alpha[index] = np.dot(xb, ij) / np.dot(ab, ij)\n            index += 1\n\n    alpha_01 = alpha[(alpha >= 0) & (alpha <= 1)]\n    alpha_avg = np.average(alpha_01)\n    alpha_std_dev = np.std(alpha_01) * 1.96\n    alpha_hist = np.histogram(alpha, 20)\n\n    return alpha, alpha_avg, alpha_std_dev, alpha_hist, alpha_01.size\n\n\n\ndef plot_fit(x, y, alpha, title, xlabel, ylabel):\n    fig, ax = plt.subplots(figsize = (plot_width_input.value, plot_height_input.value), layout = 'constrained')\n\n    ax.set_title(title, fontsize = plot_title_size_input.value)\n\n    ax.set_xlabel(xlabel, fontsize = plot_label_size_input.value)\n    ax.set_ylabel(ylabel, fontsize = plot_label_size_input.value)\n\n    ax.plot(x, y, '.')\n    ax.plot(x, alpha * x)\n\n    return fig\n\n\n\ndef plot_histogram(histogram, title, xlabel, ylabel):\n    counts = histogram[0]\n    edges = histogram[1]\n\n    fig, ax = plt.subplots(figsize = (plot_width_input.value, plot_height_input.value), layout = 'constrained')\n\n    ax.set_title(title, fontsize = plot_title_size_input.value)\n\n    ax.set_xlabel(xlabel, fontsize = plot_label_size_input.value)\n    ax.set_ylabel(ylabel, fontsize = plot_label_size_input.value)\n\n    ax.bar(edges[:-1], counts, width = np.diff(edges), edgecolor = 'black', align = 'edge')\n\n    return fig\n\n\n\ndef plot_results(event):\n    global f4ab_prime, f4xb_prime, alpha, f4ab_std, f4xb_std, alpha_std, alpha_ratio_hist\n\n    if all([el is not None for el in [f4ab_prime, f4xb_prime, alpha, f4ab_std, f4xb_std, alpha_std, alpha_ratio_hist]]):\n        f4prime_fit_fig = plot_fit(f4ab_prime, f4xb_prime, alpha, f'{hybrid} = alpha {parent1} + (1 - alpha) {parent2}', f"f4'({parent1}, {parent2}; i, j)", f"f4'({hybrid}, {parent2}; i, j)")\n        f4_fit_fig = plot_fit(f4ab_std, f4xb_std, alpha_std, f'{hybrid} = alpha {parent1} + (1 - alpha) {parent2}', f"f4({parent1}, {parent2}; i, j)", f"f4({hybrid}, {parent2}; i, j)")\n        f4_ratio_histogram_fig = plot_histogram(alpha_ratio_hist, f'{hybrid} = alpha {parent1} + (1 - alpha) {parent2}', 'f4 ratio', 'Counts')\n\n        f4prime_fit_pane.object = f4prime_fit_fig\n        f4_fit_pane.object = f4_fit_fig\n        f4_ratio_histogram_pane.object = f4_ratio_histogram_fig\n\n\n\nplot_width_input.param.watch(plot_results, 'value')\nplot_height_input.param.watch(plot_results, 'value')\nplot_title_size_input.param.watch(plot_results, 'value')\nplot_label_size_input.param.watch(plot_results, 'value')\n\n\n\ndef compute_results(event):\n    global alpha_pre_jl, cosine_pre_jl, angle_pre_jl, percentage_pre_jl, f3_test\n    global f4ab_prime, f4xb_prime, alpha, alpha_error\n    global f4ab_std, f4xb_std, alpha_std, alpha_std_error\n    global cosine_post_jl, angle_post_jl, percentage_post_jl\n    global alpha_ratio, alpha_ratio_avg, alpha_ratio_std_dev, alpha_ratio_hist, num_cases\n\n    hybrid_freqs = allele_frequencies[hybrid]\n    parent1_freqs = allele_frequencies[parent1]\n    parent2_freqs = allele_frequencies[parent2]\n    aux_freqs = [allele_frequencies[pop] for pop in auxiliaries]\n\n    alpha_pre_jl = mixing_coefficient_pre_jl(hybrid_freqs, parent1_freqs, parent2_freqs)\n\n    cosine_pre_jl, angle_pre_jl, percentage_pre_jl = admixture_angle_pre_jl(hybrid_freqs, parent1_freqs, parent2_freqs)\n\n    f3_test = f3(hybrid_freqs, parent1_freqs, parent2_freqs)\n\n    f4ab_prime, f4xb_prime = f4_prime(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs)\n    alpha, alpha_error = least_squares(f4ab_prime, f4xb_prime)\n\n    f4ab_std, f4xb_std = f4_std(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs)\n    alpha_std, alpha_std_error = least_squares(f4ab_std, f4xb_std)\n\n    cosine_post_jl, angle_post_jl, percentage_post_jl = admixture_angle_post_jl(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs)\n\n    alpha_ratio, alpha_ratio_avg, alpha_ratio_std_dev, alpha_ratio_hist, num_cases = f4_ratio(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs)\n\n    # Output\n\n    num_aux_pops = len(auxiliaries)\n    num_aux_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)\n\n    text = '### Admixture model\\n'\n    text += f'\`{hybrid} = {parent1} + {parent2}\`\\n'\n    text += f'SNPs employed: {num_valid_alleles} / {num_alleles}\\n'\n    text += f'Auxiliary populations: {num_aux_pops}\\n'\n    text += f'Auxiliary pairs: {num_aux_pairs}\\n'\n    text += '### Admixture angle\\n'\n    text += '| Calculation | Cosine | Angle (deg) | Percentage of 180 deg |\\n'\n    text += '| --- | --- | --- | --- |\\n'\n    text += f'| Pre-JL | {cosine_pre_jl:7.4f} | {angle_pre_jl:7.2f} | {percentage_pre_jl:.1%} |\\n'\n    text += f'| Post-JL | {cosine_post_jl:7.4f} | {angle_post_jl:7.2f} | {percentage_post_jl:.1%} |\\n'\n    text += '### Mixing coefficient\\n'\n    text += '| Calculation | Alpha | Error (95% CI) |\\n'\n    text += '| --- | --- | --- |\\n'\n    text += f'| Pre-JL | {alpha_pre_jl:6.4f} | - |\\n'\n    text += f'| Post-JL (f4-prime, renormalized) | {alpha:6.4f} | {alpha_error:6.4f} |\\n'\n    text += f'| Post-JL NR (f4, standard) | {alpha_std:6.4f} | {alpha_std_error:6.4f} |\\n'\n    text += '### f4-ratio and f3 test\\n'\n    text += f'f4-ratio average if [0, 1]: {alpha_ratio_avg:6.4f} +/- {alpha_ratio_std_dev:6.4f} (95% CI), {num_cases} cases\\n'\n    text += f'Standard admixture test: f3(c1, c2; x) < 0 ? {f3_test:8.6f}'\n\n    resulting_data_output.object = text\n\n    plot_results(None)\n\n    results_download_button.disabled = False\n    f4_points_download_button.disabled = False\n\n\n\npn.bind(compute_results, compute_results_button, watch = True)\n\n\n\ndef save_population_allele_frequencies():\n    file = StringIO()\n\n    pops_width = max([len(name) for name in sel_pops])\n    prec = 6\n    col_width = max(prec + 7, pops_width)\n\n    headers_format = ' '.join([f'{{{i}:^{col_width}}}' for i, pop in enumerate(sel_pops)])\n    headers = headers_format.format(*sel_pops)\n    file.write(headers + '\\n')\n\n    row_format = ' '.join([f'{{{i}: {col_width}.{prec}E}}' for i, pop in enumerate(sel_pops)])\n\n    for allele_index in range(num_valid_alleles):\n        row = [freqs[allele_index] for pop, freqs in allele_frequencies.items()]\n        file.write(row_format.format(*row) + '\\n')\n\n    file.seek(0)\n\n    now = datetime.now()\n    name = now.strftime("frequencies_%Y-%m-%d_%Hh%Mm%Ss")\n    freqs_download_button.filename = name + '.dat'\n\n    return file\n\n\n\nfreqs_download_button.callback = save_population_allele_frequencies\n\n\n\ndef save_f4_points():\n    file = StringIO()\n\n    aux_pops_width = max([len(name) for name in auxiliaries])\n    prec = 6\n    col_width = prec + 7\n\n    headers = '{0:^{col_width}} {1:^{col_width}} {2:^{col_width}} {3:^{col_width}} {4:^{col_width}} {5:^{aux_pops_width}} {6:^{aux_pops_width}}'.format('f4primeAB', 'f4primeXB', 'f4AB', 'f4XB', 'f4-ratio', 'Aux1', 'Aux2', col_width = col_width, aux_pops_width = aux_pops_width)\n    file.write(headers + '\\n')\n\n    num_aux_pops = len(sel_pops[3:])\n\n    index = 0\n\n    for i in range(num_aux_pops):\n        for j in range(i + 1, num_aux_pops):\n            row = '{0: {col_width}.{prec}E} {1: {col_width}.{prec}E} {2: {col_width}.{prec}E} {3: {col_width}.{prec}E} {4: {col_width}.{prec}E} {5:{aux_pops_width}} {6:{aux_pops_width}}'.format(f4ab_prime[index], f4xb_prime[index], f4ab_std[index], f4xb_std[index], alpha_ratio[index], auxiliaries[i], auxiliaries[j], prec = prec, col_width = col_width, aux_pops_width = aux_pops_width)\n            file.write(row + '\\n')\n            index += 1\n\n    file.seek(0)\n\n    now = datetime.now()\n    name = now.strftime("f4_points_%Y-%m-%d_%Hh%Mm%Ss")\n    f4_points_download_button.filename = name + '.dat'\n\n    return file\n\n\n\nf4_points_download_button.callback = save_f4_points\n\n\n\ndef save_admixture_data():\n    num_aux_pops = len(auxiliaries)\n    num_aux_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)\n\n    file = StringIO()\n\n    file.write(f'Admixture model: {hybrid} = {parent1} + {parent2}\\n')\n    file.write(f'SNPs employed: {num_valid_alleles} / {num_alleles}\\n')\n    file.write(f'Auxiliary populations: {num_aux_pops}\\n')\n    file.write(f'Auxiliary pairs: {num_aux_pairs}\\n')\n    file.write(f'Cos pre-JL:  {cosine_pre_jl:7.4f} ---> Angle pre-JL:  {angle_pre_jl:7.2f} deg vs 180 deg: {percentage_pre_jl:.1%}\\n')\n    file.write(f'Cos post-JL: {cosine_post_jl:7.4f} ---> Angle post-JL: {angle_post_jl:7.2f} deg vs 180 deg: {percentage_post_jl:.1%}\\n')\n    file.write(f'Alpha pre-JL:     {alpha_pre_jl:6.4f}\\n')\n    file.write(f'Alpha post-JL:    {alpha:6.4f} +/- {alpha_error:6.4f} (95% CI) (f4-prime, renormalized)\\n')\n    file.write(f'Alpha NR post-JL: {alpha_std:6.4f} +/- {alpha_std_error:6.4f} (95% CI) (f4, standard)\\n')\n    file.write(f'f4-ratio average if [0, 1]: {alpha_ratio_avg:6.4f} +/- {alpha_ratio_std_dev:6.4f} (95% CI), {num_cases} cases\\n')\n    file.write(f'Standard admixture test: f3(c1, c2; x) < 0 ? {f3_test:8.6f}')\n\n    file.seek(0)\n\n    now = datetime.now()\n    name = now.strftime("results_%Y-%m-%d_%Hh%Mm%Ss")\n    results_download_button.filename = name + '.dat'\n\n    return file\n\n\n\nresults_download_button.callback = save_admixture_data\n\n\n\ndef admixture_model_card_width():\n    return hybrid_select.width + 2 * hybrid_select.margin[1] + parent1_select.width + 2 * parent1_select.margin[1] + parent2_select.width + 2 * parent2_select.margin[1] + aux_pops_select.width + 2 * aux_pops_select.margin[1]\n\ndef admixture_model_card_height():\n    print(hybrid_select.height)\n    return hybrid_select.height + 2 * hybrid_select.margin[0]\n\n\n\n# Cards and their layouts\n\ncard_margin = 10\n\navail_pops_box = pn.WidgetBox('### Available populations', avail_pops_filter, avail_pops_table, margin = card_margin)\nsel_pops_box = pn.WidgetBox('### Selected populations', reset_sel_pops_button, sel_pops_table, margin = card_margin)\ncomp_box = pn.WidgetBox('### Computation', num_procs_input, compute_freqs_button, freqs_download_button, margin = card_margin)\n\nadmixture_model_layout = pn.FlexBox(hybrid_select, parent1_select, parent2_select, aux_pops_select, flex_direction = 'row', justify_content = 'space-evenly')\nadmixture_model_card = pn.Card(admixture_model_layout, title = 'Admixture model', collapsible = False, margin = card_margin, styles = {'width': 'fit-content'})\n\nactions_card = pn.Card(compute_results_button, f4_points_download_button, results_download_button, title = 'Actions', collapsible = False, margin = card_margin, styles = {'width': 'fit-content'})\nplot_options_card = pn.Card(plot_width_input, plot_height_input, plot_title_size_input, plot_label_size_input, title = 'Plot options', collapsible = False, margin = card_margin, styles = {'width': 'fit-content'})\nresulting_data_card = pn.Card(resulting_data_output, title = 'Results', collapsible = False, margin = card_margin, styles = {'width': 'fit-content'})\n\nflex_box1 = pn.FlexBox(actions_card, plot_options_card, flex_direction = 'row', justify_content = 'space-evenly', margin = card_margin, styles = {'width': 'fit-content'})\n\ncol = pn.Column(admixture_model_card, flex_box1, styles = {'width': 'fit-content'})\n\nflex_box2 = pn.FlexBox(col, resulting_data_card, flex_direction = 'row', justify_content = 'space-evenly', margin = card_margin, styles = {'width': 'fit-content'})\n\nplots_layout = pn.FlexBox(f4prime_fit_pane, f4_fit_pane, f4_ratio_histogram_pane, flex_direction = 'row', justify_content = 'space-evenly')\nplots_card = pn.Card(plots_layout, title = 'Plots', collapsible = False, margin = card_margin, styles = {'width': 'fit-content'})\n\n# Layouts\n\ninput_files_layout = pn.Column(file_selector, load_files_button, name = 'Inputs', sizing_mode = 'stretch_width')\nselect_pops_layout = pn.FlexBox(avail_pops_box, sel_pops_box, comp_box, name = 'Select populations', flex_direction = 'row')\nresults_layout = pn.Column(flex_box2, plots_card, name = 'Results')\n\n# Tabs\n\ntabs = pn.Tabs(input_files_layout, select_pops_layout, results_layout)\n\n# Main layout\n\nmain_layout = pn.Column(alert_pane, tabs, scroll = True)\n\n# Template\n\ntemplate = pn.template.VanillaTemplate(title = 'Mixtum: The geometry of admixture in population genetics')\ntemplate.main.append(main_layout)\ntemplate.servable()\n\n\nawait write_doc()
  `

  try {
    const [docs_json, render_items, root_ids] = await self.pyodide.runPythonAsync(code)
    self.postMessage({
      type: 'render',
      docs_json: docs_json,
      render_items: render_items,
      root_ids: root_ids
    })
  } catch(e) {
    const traceback = `${e}`
    const tblines = traceback.split('\n')
    self.postMessage({
      type: 'status',
      msg: tblines[tblines.length-2]
    });
    throw e
  }
}

self.onmessage = async (event) => {
  const msg = event.data
  if (msg.type === 'rendered') {
    self.pyodide.runPythonAsync(`
    from panel.io.state import state
    from panel.io.pyodide import _link_docs_worker

    _link_docs_worker(state.curdoc, sendPatch, setter='js')
    `)
  } else if (msg.type === 'patch') {
    self.pyodide.globals.set('patch', msg.patch)
    self.pyodide.runPythonAsync(`
    from panel.io.pyodide import _convert_json_patch
    state.curdoc.apply_json_patch(_convert_json_patch(patch), setter='js')
    `)
    self.postMessage({type: 'idle'})
  } else if (msg.type === 'location') {
    self.pyodide.globals.set('location', msg.location)
    self.pyodide.runPythonAsync(`
    import json
    from panel.io.state import state
    from panel.util import edit_readonly
    if state.location:
        loc_data = json.loads(location)
        with edit_readonly(state.location):
            state.location.param.update({
                k: v for k, v in loc_data.items() if k in state.location.param
            })
    `)
  }
}

startApplication()
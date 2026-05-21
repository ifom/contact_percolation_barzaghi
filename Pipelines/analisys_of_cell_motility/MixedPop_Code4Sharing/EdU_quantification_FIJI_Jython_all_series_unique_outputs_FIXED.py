#@ File    (label="Cartella input: contiene le 3 cartelle/replicati", style="directory") input_folder
#@ File    (label="Cartella output risultati", style="directory") output_folder
#@ String  (label="Estensione file", value=".nd") ext
#@ Integer (label="Canale DAPI / nuclei", value=1) CH_DAPI
#@ Integer (label="Canale EdU", value=3) CH_EDU
#@ String  (label="Threshold mode: manual oppure negative_percentile", value="manual") threshold_mode
#@ Float   (label="Manual EdU threshold su EdU_Mean", value=500.0) manual_threshold
#@ String  (label="Substring controllo negativo nel nome condizione", value="negative") negative_condition_text
#@ Float   (label="Percentile controllo negativo per soglia", value=99.0) negative_percentile
#@ Boolean (label="Analizza tutte le serie/FOV nei file Bio-Formats", value=True) process_all_series
#@ Integer (label="Max file da processare; 0 = tutti", value=0) max_files_total
#@ Integer (label="Max serie/FOV per file; 0 = tutte", value=0) max_series_per_file
#@ Integer (label="Resize StarDist width", value=256) stardist_resize_width
#@ Integer (label="Resize StarDist height", value=256) stardist_resize_height
#@ Float   (label="Subtract background EdU rolling ball; 0 = no", value=0.0) edu_background_rolling
#@ Float   (label="Area nuclei minima; 0 = nessun filtro", value=0.0) min_nucleus_area
#@ Float   (label="Area nuclei massima; 0 = nessun filtro", value=0.0) max_nucleus_area
#@ Boolean (label="Salva immagini annotate", value=True) save_annotated_images
#@ Boolean (label="Salva ROI zip", value=True) save_roi_zip
#@ CommandService command

# EdU_quantification_FIJI_Jython_all_series_unique_outputs_FIXED.py
#
# Quantificazione EdU cellula-per-cellula da file Bio-Formats (.nd).
# Modifica principale rispetto alla versione precedente:
#   - immagini annotate e ROI zip hanno nomi univoci tra replicati/FOV.
#   - fix della funzione save_annotation_and_rois(): non usa piu' variabili globali inesistenti.

from ij import IJ
from ij.measure import ResultsTable, Measurements
from ij.plugin import Duplicator, RoiScaler
from ij.plugin.frame import RoiManager
from de.csbdresden.stardist import StarDist2D
from loci.plugins import BF
import os
import math

# ImporterOptions vive in un package chiamato "in", che e' una keyword Python.
ImporterOptions = __import__('loci.plugins.in', globals(), locals(), ['ImporterOptions']).ImporterOptions


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def safe_name(text):
    bad = ['/', '\\', ':', '*', '?', '"', '<', '>', '|', ' ', '\t', '\n', '\r']
    out = str(text)
    for b in bad:
        out = out.replace(b, '_')
    while '__' in out:
        out = out.replace('__', '_')
    out = out.strip('_')
    if out == '':
        out = 'unnamed'
    return out


def basename_without_ext(filename):
    base = os.path.basename(filename)
    if base.lower().endswith(ext.lower()):
        return base[:-len(ext)]
    return os.path.splitext(base)[0]


def infer_replicate(src_dir, current_dir):
    rel = os.path.relpath(current_dir, src_dir)
    if rel == '.' or rel == os.curdir:
        return 'Replicate_Unknown'
    parts = rel.split(os.sep)
    if len(parts) > 0 and parts[0] != '':
        return parts[0]
    return 'Replicate_Unknown'


def make_unique_output_base(src_dir, current_dir, filename, series_index):
    replicate = infer_replicate(src_dir, current_dir)
    rel_path = os.path.relpath(current_dir, src_dir)
    if rel_path == '.' or rel_path == os.curdir:
        rel_path = 'root'

    rel_path_safe = safe_name(rel_path.replace(os.sep, '_'))
    condition_safe = safe_name(basename_without_ext(filename))
    replicate_safe = safe_name(replicate)

    # Nome univoco tra replicati, sottocartelle, file e FOV/serie.
    # Esempio:
    # 2026_04_15__2026_04_15__RAB5A_KO_STAT2A__series_14
    return safe_name(
        replicate_safe +
        '__' + rel_path_safe +
        '__' + condition_safe +
        '__series_' + str(int(series_index))
    )


def get_roi_manager(reset):
    rm = RoiManager.getInstance()
    if rm is None:
        rm = RoiManager()
    if reset:
        rm.reset()
    return rm


def close_if_not_none(imp):
    try:
        if imp is not None:
            imp.close()
    except:
        pass


def percentile(values, pct):
    vals = [float(v) for v in values]
    vals.sort()
    n = len(vals)
    if n == 0:
        return None
    if n == 1:
        return vals[0]
    if pct <= 0:
        return vals[0]
    if pct >= 100:
        return vals[-1]
    pos = (pct / 100.0) * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac


def mean(values):
    if len(values) == 0:
        return 0.0
    return sum(values) / float(len(values))


def median(values):
    if len(values) == 0:
        return 0.0
    vals = [float(v) for v in values]
    vals.sort()
    n = len(vals)
    mid = n // 2
    if n % 2 == 1:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2.0


def stdev(values):
    if len(values) < 2:
        return 0.0
    m = mean(values)
    ss = 0.0
    for v in values:
        ss += (v - m) * (v - m)
    return math.sqrt(ss / float(len(values) - 1))


def sem(values):
    if len(values) < 2:
        return 0.0
    return stdev(values) / math.sqrt(float(len(values)))


def open_bioformats_imageplus_list(path):
    options = ImporterOptions()
    options.setId(path)
    options.setAutoscale(True)
    options.setOpenAllSeries(process_all_series)
    try:
        options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
    except:
        pass
    imps = BF.openImagePlus(options)
    return imps


def segment_nuclei_stardist(imp):
    # Segmenta sul canale DAPI, eventualmente ridimensionato per velocizzare StarDist.
    # Le ROI generate sul resize vengono riportate alle dimensioni originali con scale_x/scale_y.
    rm = get_roi_manager(True)
    dupl = Duplicator()

    dup_nuclei = dupl.run(imp, CH_DAPI, CH_DAPI, 1, 1, 1, 1)
    orig_width = dup_nuclei.getWidth()
    orig_height = dup_nuclei.getHeight()

    scale_x = 1.0
    scale_y = 1.0

    if stardist_resize_width > 0 and stardist_resize_height > 0:
        ip_nuclei = dup_nuclei.getProcessor().resize(stardist_resize_width, stardist_resize_height)
        dup_nuclei.setProcessor(ip_nuclei)
        scale_x = float(orig_width) / float(stardist_resize_width)
        scale_y = float(orig_height) / float(stardist_resize_height)

    label = None
    try:
        res = command.run(StarDist2D, False,
            'input', dup_nuclei,
            'modelChoice', 'Versatile (fluorescent nuclei)',
            'normalizeInput', True,
            'outputType', 'Both'
        ).get()
        label = res.getOutput('label')
    except Exception as e:
        print('ERRORE StarDist su {}: {}'.format(imp.getTitle(), e))

    close_if_not_none(dup_nuclei)

    rm = get_roi_manager(False)
    n = rm.getCount()

    # Riporta ROI alla dimensione originale se e' stato fatto resize.
    if n > 0 and (abs(scale_x - 1.0) > 0.0001 or abs(scale_y - 1.0) > 0.0001):
        for r in range(n):
            roi_n = rm.getRoi(r)
            roi_scaled = RoiScaler.scale(roi_n, scale_x, scale_y, 0)
            rm.setRoi(roi_scaled, r)

    close_if_not_none(label)
    return rm.getCount()


def measure_cells_in_imp(imp, src_dir, current_dir, filename, series_index, series_title):
    rows = []
    replicate = infer_replicate(src_dir, current_dir)
    condition = basename_without_ext(filename)

    dupl = Duplicator()
    rm = get_roi_manager(True)
    nuclei_count_raw = segment_nuclei_stardist(imp)
    rm = get_roi_manager(False)

    if nuclei_count_raw == 0:
        print('ATTENZIONE: nessun nucleo trovato in file={}, serie={}'.format(filename, series_index))
        return rows

    edu_imp = dupl.run(imp, CH_EDU, CH_EDU, 1, 1, 1, 1)
    dapi_imp = dupl.run(imp, CH_DAPI, CH_DAPI, 1, 1, 1, 1)

    if edu_background_rolling > 0:
        try:
            IJ.run(edu_imp, 'Subtract Background...', 'rolling={}'.format(edu_background_rolling))
        except Exception as e:
            print('ATTENZIONE: background subtraction EdU fallita: {}'.format(e))

    meas = Measurements.AREA + Measurements.MEAN + Measurements.MIN_MAX + Measurements.MEDIAN

    valid_cell_id = 0
    for i in range(rm.getCount()):
        roi = rm.getRoi(i)

        edu_imp.setRoi(roi)
        dapi_imp.setRoi(roi)

        edu_stats = edu_imp.getStatistics(meas)
        dapi_stats = dapi_imp.getStatistics(meas)

        area = float(edu_stats.area)
        if min_nucleus_area > 0 and area < min_nucleus_area:
            continue
        if max_nucleus_area > 0 and area > max_nucleus_area:
            continue

        valid_cell_id += 1
        edu_mean = float(edu_stats.mean)
        edu_median = float(edu_stats.median)
        edu_intden = area * edu_mean
        dapi_mean = float(dapi_stats.mean)
        dapi_median = float(dapi_stats.median)
        dapi_intden = area * dapi_mean

        row = {
            'Replicate': replicate,
            'Condition': condition,
            'Filename': filename,
            'Relative_Path': os.path.relpath(current_dir, src_dir),
            'Series_Index': int(series_index),
            'Series_Title': series_title,
            'Cell_ID': int(valid_cell_id),
            'Nuclei_Count_Raw_FOV': int(nuclei_count_raw),
            'Nuclei_Count_Valid_FOV': 0,
            'Nucleus_Area': area,
            'EdU_Mean': edu_mean,
            'EdU_Median': edu_median,
            'EdU_IntDen': edu_intden,
            'DAPI_Mean': dapi_mean,
            'DAPI_Median': dapi_median,
            'DAPI_IntDen': dapi_intden,
            'EdU_Positive': 0,
            'EdU_Threshold_Used': 0.0
        }
        rows.append(row)

    for row in rows:
        row['Nuclei_Count_Valid_FOV'] = int(valid_cell_id)

    close_if_not_none(edu_imp)
    close_if_not_none(dapi_imp)
    return rows


def save_annotation_and_rois(imp, dst_dir, src_dir, current_dir, filename, series_index):
    rm = RoiManager.getInstance()
    if rm is None or rm.getCount() == 0:
        return

    base = make_unique_output_base(src_dir, current_dir, filename, series_index)

    if save_roi_zip:
        roi_dir = os.path.join(dst_dir, 'roi_zip')
        ensure_dir(roi_dir)
        roi_path = os.path.join(roi_dir, base + '_nuclei_rois.zip')
        try:
            rm.runCommand('Save', roi_path)
        except Exception as e:
            print('ATTENZIONE: impossibile salvare ROI zip {}: {}'.format(roi_path, e))

    if save_annotated_images:
        ann_dir = os.path.join(dst_dir, 'annotated_images')
        ensure_dir(ann_dir)
        ann_path = os.path.join(ann_dir, base + '_annotated.tif')
        try:
            imp.show()
            try:
                imp.setDisplayMode(IJ.COLOR)
            except:
                pass
            rm.runCommand(imp, 'Show All with labels')
            IJ.saveAs(imp, 'Tiff', ann_path)
            imp.hide()
        except Exception as e:
            print('ATTENZIONE: impossibile salvare immagine annotata {}: {}'.format(ann_path, e))


def process_file(src_dir, dst_dir, current_dir, filename):
    all_rows = []
    path = os.path.join(current_dir, filename)
    print('\nApro file: {}'.format(path))

    imps = None
    try:
        imps = open_bioformats_imageplus_list(path)
    except Exception as e:
        print('ERRORE apertura Bio-Formats file {}: {}'.format(path, e))
        return all_rows

    if imps is None or len(imps) == 0:
        print('ATTENZIONE: nessuna serie aperta per {}'.format(path))
        return all_rows

    n_series = len(imps)
    if max_series_per_file > 0 and max_series_per_file < n_series:
        n_series = max_series_per_file

    for s in range(n_series):
        imp = imps[s]
        if imp is None:
            continue
        series_index = s + 1
        series_title = imp.getTitle()
        print('  Serie/FOV {} di {}: {}'.format(series_index, n_series, series_title))

        try:
            rows = measure_cells_in_imp(imp, src_dir, current_dir, filename, series_index, series_title)
            all_rows.extend(rows)
            save_annotation_and_rois(imp, dst_dir, src_dir, current_dir, filename, series_index)
        except Exception as e:
            print('ERRORE processing file={}, serie={}: {}'.format(filename, series_index, e))

        close_if_not_none(imp)
        rm = RoiManager.getInstance()
        if rm is not None:
            rm.reset()

    # Chiude eventuali serie non processate quando max_series_per_file limita l'analisi.
    try:
        for s in range(n_series, len(imps)):
            close_if_not_none(imps[s])
    except:
        pass

    return all_rows


def collect_all_measurements():
    src_dir = input_folder.getAbsolutePath()
    dst_dir = output_folder.getAbsolutePath()
    ensure_dir(dst_dir)

    all_rows = []
    file_counter = 0

    for root, directories, filenames in os.walk(src_dir):
        directories.sort()
        filenames.sort()
        for filename in filenames:
            if not filename.lower().endswith(ext.lower()):
                continue
            if max_files_total > 0 and file_counter >= max_files_total:
                return all_rows
            file_counter += 1
            rows = process_file(src_dir, dst_dir, root, filename)
            all_rows.extend(rows)

    return all_rows


def determine_threshold(all_rows):
    mode = str(threshold_mode).strip().lower()
    if mode.startswith('negative'):
        neg_text = str(negative_condition_text).strip().lower()
        neg_values = []
        for row in all_rows:
            cond = str(row['Condition']).lower()
            if neg_text != '' and cond.find(neg_text) >= 0:
                neg_values.append(float(row['EdU_Mean']))
        thr = percentile(neg_values, float(negative_percentile))
        if thr is None:
            print('ATTENZIONE: nessuna cellula trovata per controllo negativo contenente "{}". Uso soglia manuale.'.format(negative_condition_text))
            return float(manual_threshold), 'manual_fallback'
        print('Soglia EdU calcolata dal controllo negativo: percentile {} = {}'.format(negative_percentile, thr))
        return float(thr), 'negative_percentile'

    return float(manual_threshold), 'manual'


def apply_threshold(all_rows, threshold):
    for row in all_rows:
        row['EdU_Threshold_Used'] = float(threshold)
        if float(row['EdU_Mean']) >= float(threshold):
            row['EdU_Positive'] = 1
        else:
            row['EdU_Positive'] = 0


def add_row_to_rt(rt, row, columns):
    rt.incrementCounter()
    for c in columns:
        rt.addValue(c, row[c])


def save_cell_table(all_rows, dst_dir):
    columns = [
        'Replicate', 'Condition', 'Filename', 'Relative_Path', 'Series_Index', 'Series_Title',
        'Cell_ID', 'Nuclei_Count_Raw_FOV', 'Nuclei_Count_Valid_FOV', 'Nucleus_Area',
        'EdU_Mean', 'EdU_Median', 'EdU_IntDen',
        'DAPI_Mean', 'DAPI_Median', 'DAPI_IntDen',
        'EdU_Threshold_Used', 'EdU_Positive'
    ]
    rt = ResultsTable()
    for row in all_rows:
        add_row_to_rt(rt, row, columns)
    rt.show('EdU_cells_results')
    out = os.path.join(dst_dir, 'EdU_cells_results.csv')
    rt.save(out)
    print('Salvato: {}'.format(out))


def summarize_group(rows, group_keys):
    groups = {}
    for row in rows:
        key = []
        for k in group_keys:
            key.append(row[k])
        key = tuple(key)
        if key not in groups:
            groups[key] = []
        groups[key].append(row)

    out_rows = []
    for key in sorted(groups.keys()):
        gr = groups[key]
        edu_mean_values = [float(r['EdU_Mean']) for r in gr]
        edu_median_values = [float(r['EdU_Median']) for r in gr]
        edu_intden_values = [float(r['EdU_IntDen']) for r in gr]
        area_values = [float(r['Nucleus_Area']) for r in gr]
        pos_values = [float(r['EdU_Positive']) for r in gr]

        out = {}
        for i in range(len(group_keys)):
            out[group_keys[i]] = key[i]

        n_cells = len(gr)
        n_pos = int(sum(pos_values))
        pct_pos = 0.0
        if n_cells > 0:
            pct_pos = 100.0 * float(n_pos) / float(n_cells)

        out['N_Cells'] = int(n_cells)
        out['N_EdU_Positive'] = int(n_pos)
        out['Percent_EdU_Positive'] = float(pct_pos)
        out['EdU_Mean_Mean'] = mean(edu_mean_values)
        out['EdU_Mean_Median'] = median(edu_mean_values)
        out['EdU_Mean_SD'] = stdev(edu_mean_values)
        out['EdU_Mean_SEM'] = sem(edu_mean_values)
        out['EdU_Median_Mean'] = mean(edu_median_values)
        out['EdU_IntDen_Mean'] = mean(edu_intden_values)
        out['Nucleus_Area_Mean'] = mean(area_values)
        out['EdU_Threshold_Used'] = float(gr[0]['EdU_Threshold_Used'])
        out_rows.append(out)

    return out_rows


def save_summary_table(rows, dst_dir, filename, group_keys):
    summary_rows = summarize_group(rows, group_keys)
    columns = []
    for k in group_keys:
        columns.append(k)
    columns.extend([
        'N_Cells', 'N_EdU_Positive', 'Percent_EdU_Positive',
        'EdU_Mean_Mean', 'EdU_Mean_Median', 'EdU_Mean_SD', 'EdU_Mean_SEM',
        'EdU_Median_Mean', 'EdU_IntDen_Mean', 'Nucleus_Area_Mean', 'EdU_Threshold_Used'
    ])

    rt = ResultsTable()
    for row in summary_rows:
        add_row_to_rt(rt, row, columns)
    rt.show(filename.replace('.csv', ''))
    out = os.path.join(dst_dir, filename)
    rt.save(out)
    print('Salvato: {}'.format(out))


def save_threshold_info(dst_dir, threshold, threshold_source):
    rt = ResultsTable()
    rt.incrementCounter()
    rt.addValue('Threshold_Mode_Requested', str(threshold_mode))
    rt.addValue('Threshold_Source_Used', str(threshold_source))
    rt.addValue('EdU_Threshold_Used', float(threshold))
    rt.addValue('Manual_Threshold', float(manual_threshold))
    rt.addValue('Negative_Condition_Text', str(negative_condition_text))
    rt.addValue('Negative_Percentile', float(negative_percentile))
    out = os.path.join(dst_dir, 'EdU_threshold_info.csv')
    rt.save(out)
    print('Salvato: {}'.format(out))


def main():
    dst_dir = output_folder.getAbsolutePath()
    ensure_dir(dst_dir)

    rm = RoiManager.getInstance()
    if rm is not None:
        rm.close()
    RoiManager()

    print('\n=== Inizio analisi EdU ===')
    print('Input: {}'.format(input_folder.getAbsolutePath()))
    print('Output: {}'.format(dst_dir))
    print('Estensione: {}'.format(ext))
    print('DAPI channel: {}'.format(CH_DAPI))
    print('EdU channel: {}'.format(CH_EDU))

    all_rows = collect_all_measurements()

    if len(all_rows) == 0:
        print('Nessuna cellula misurata. Controlla input, estensione, canali e StarDist.')
        return

    threshold, threshold_source = determine_threshold(all_rows)
    apply_threshold(all_rows, threshold)

    save_cell_table(all_rows, dst_dir)

    save_summary_table(
        all_rows,
        dst_dir,
        'EdU_summary_by_FOV.csv',
        ['Replicate', 'Condition', 'Filename', 'Series_Index', 'Series_Title']
    )

    save_summary_table(
        all_rows,
        dst_dir,
        'EdU_summary_by_replicate_condition.csv',
        ['Replicate', 'Condition']
    )

    save_summary_table(
        all_rows,
        dst_dir,
        'EdU_summary_by_condition_preview.csv',
        ['Condition']
    )

    save_threshold_info(dst_dir, threshold, threshold_source)

    rm = RoiManager.getInstance()
    if rm is not None:
        rm.close()

    print('\n=== Analisi completata ===')
    print('Cellule misurate: {}'.format(len(all_rows)))
    print('Soglia EdU usata: {}'.format(threshold))
    print('Risultati principali:')
    print('  - EdU_cells_results.csv')
    print('  - EdU_summary_by_FOV.csv')
    print('  - EdU_summary_by_replicate_condition.csv')
    print('  - EdU_summary_by_condition_preview.csv')


main()

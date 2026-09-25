import json
import pytest
from pn_fusion.validation import check_data_integrity,check_models,check_fingerprints,check_trajectories,audit_archived_work
from pn_fusion.reporting import export_tables
from pn_fusion.io import ROOT,write_csv
from pn_fusion.benchmarks import groups,run_benchmark

def test_archive_numeric_integrity():
    assert check_data_integrity()['unaltered_numeric_cells']==23302

def test_archived_model_values():
    assert check_models()['scalar_checks']==1682

def test_canonical_matrix_fingerprints():
    assert check_fingerprints()['matrix_fingerprints']==15

def test_all_archived_trajectories():
    assert check_trajectories()['archived_configurations_checked']==606

def test_known_work_discrepancy_is_not_hidden(tmp_path):
    report=audit_archived_work(tmp_path)
    assert report['discrepant_values']==320
    assert report['affected_datasets']==['fusion_end_to_end_results.csv','fusion_stationary_results.csv']
    assert report['archived_data_modified'] is False

def test_all_manuscript_table_exports(tmp_path):
    result=export_tables(tmp_path)
    assert len(result)==22
    assert len(list(tmp_path.glob('table_*.csv')))==22

def test_reference_outputs_are_protected():
    with pytest.raises(ValueError):write_csv(ROOT/'data/reference/SHOULD_NOT_EXIST.csv',[{'x':1}])
    assert not (ROOT/'data/reference/SHOULD_NOT_EXIST.csv').exists()

def test_configuration_inventory():
    expected={'initialization':24,'sweep':234,'p5-m6':90,'p7-m8':90,'hammerstein':12,'fusion':144,'stationary':16,'phat2-p10':36}
    assert sum(expected.values())==646
    for suite,count in expected.items():assert sum(len(g[4]) for g in groups(suite))==count

def test_raw_sample_smoke(tmp_path):
    report=run_benchmark('phat2-p10',tmp_path,repeats=3,warmup=1,threads=1,quick=True)
    assert report['runs'][0]['raw_samples']==6
    assert (tmp_path/'phat2-p10_samples.csv').is_file()
    env=json.loads((tmp_path/'environment.json').read_text())
    assert env['quick_mode'] is True

#!/usr/bin/env python3
"""Reanalyse strictly transcribed real/published observational tables.
No synthetic scientific data, random resampling, p-values, or model refitting.
"""
from __future__ import annotations
import csv, hashlib, json, math, pathlib, sys, datetime
from dataclasses import asdict

ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from dcae import completeness_proxy, ScaleContract, admissibility


def rows(p):
    with (ROOT/'data'/p).open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def compute():
    walk=rows('cosmos_2026_walk_table2.csv')
    loss=rows('cosmos_2026_loss_table3.csv')
    comps=rows('cosmos_2026_components_table4.csv')
    boundary=rows('cosmos_2026_boundary.csv')
    scsn=rows('scsn_degradation.csv')
    assert len(walk)==15 and len(loss)==9 and len(comps)==2 and len(boundary)==5 and len(scsn)==9
    w=next(x for x in walk if x['walkers']=='400')
    l=next(x for x in loss if x['p_loss']=='0.1')
    d_w=float(w['coverage']); d_l=float(l['coverage'])
    # A coverage-only, geometry-agnostic *diagnostic*, not full DCAE.
    gap_w=1-math.sqrt(d_w)
    gap_l=1-math.sqrt(d_l)
    comparison={
        'source':'Imbrisak_Tisanić_2026_Tables_2_3_4',
        'experimental_unit':'published_summary_over_subsamples_of_single_COSMOS_2000_galaxy_patch',
        'walk_coverage':d_w,'loss_coverage':d_l,
        'coverage_difference_pp':100*(d_l-d_w),
        'geometric_gap_0_1_walk':gap_w,
        'geometric_gap_0_1_loss':gap_l,
        'geometric_gap_difference_abs':abs(gap_w-gap_l),
        'tolerance_contrasts':{},
        'component_relative_rmse':{},
        'reason_not_full_DCAE':'coverage is within a fixed field; L_F(s) across physically distinct scales is not measured',
    }
    for e in ('1.10','1.54','2.03','2.50'):
        hw=float(w['H_'+e]);sl=float(l['S_'+e]);
        comparison['tolerance_contrasts'][e]={
            'walk_success':hw, 'independent_loss_success':sl,
            'walk_minus_loss_pp':100*(hw-sl),
            'preferred_under_this_threshold': ('walk' if hw>sl else 'independent_loss' if hw<sl else 'tie'),
            'values_are_published_rounded':True}
    for compname in ('delta_omega','f_hd','L_hd'):
        a=float(comps[0]['rel_rmse_'+compname]);b=float(comps[1]['rel_rmse_'+compname])
        comparison['component_relative_rmse'][compname]={
            'walk':a,'independent_loss':b,'walk_over_loss':a/b,
            'lower_error_mechanism':'walk' if a<b else 'independent_loss'}
    # The count of excluded boundary carriers is a STATISTIC-SPECIFIC computational
    # support indicator; it is not the physical observer's angular mask.
    two=next(x for x in boundary if x['galaxies']=='2000')
    full=next(x for x in boundary if x['galaxies']=='94043')
    c0=int(two['interior_exact'])/int(two['galaxies'])
    c1=int(full['interior_exact'])/int(full['galaxies'])
    support_asym=completeness_proxy(c0,c1,d=2)
    field={
        'status':'SURROGATE_ONLY_NOT_DCAE_VALIDATION',
        'small_N':2000,'large_N':94043,
        'small_area_deg2_approx':float(two['field_area_deg2']),
        'large_area_deg2_approx':float(full['field_area_deg2']),
        'angular_linear_ratio_approx':math.sqrt(float(full['field_area_deg2'])/float(two['field_area_deg2'])),
        'interior_frac_small':c0,'interior_frac_large':c1,
        'interior_increase_percentage_points':100*(c1-c0),
        'carrier_support_asymmetry':asdict(support_asym),
        'why_not_observational_extent':'eligibility of bounded Voronoi cells is estimator support, NOT instrument sky mask or cosmic physical target coverage',
    }
    levels=[float(x['sigma_h_km']) for x in scsn]
    pbnds=[float(x['Pbnd_pct']) for x in scsn]
    scsn_check={
        'status':'SEPARATE_INSTRUMENTAL_STRESS_REANALYSIS_NOT_PHYSICAL_SCALE_SERIES',
        'points':len(scsn), 'sigma_range_km':[min(levels),max(levels)],
        'reported_posterior_boundary_mass_pct_first_last':[pbnds[0],pbnds[-1]],
        'strictly_increasing_reported_boundary_mass':all(pbnds[i]<pbnds[i+1] for i in range(len(pbnds)-1)),
        'physical_extent_measured_at_each_sigma':False,
        'does_not_test_DCAE_scale_growth':True,
        'external_prospective_contradiction':'Sysmic-X2 reported FAIL for threshold mechanism; see observational-degeneration-synergy/CURRENT_EVIDENCE.md',
        'source_git_blob_sha1':'54225872f65f7949527a6a8ae049868fd275e3dd',
    }
    # Deliberately do not fabricate L_F/L_O from sample count, geolocation
    # uncertainty, or boundary posterior mass.
    contract=ScaleContract(
        domain='cosmic_web',target_observable='filament_connected_extent',
        coordinate_frame='comoving_3D_not_yet_specified_for_this_table',length_unit='cMpc',
        physical_extent_definition='target_connected_extent_or_morphological_scale_required',
        observational_extent_definition='matched_selection_kernel_and_radial_anglar_resolution_required',
        source_identifier='arxiv:2609.26140 Table 2-4 reported aggregate tables',
        source_status='SECONDARY_PUBLISHED_SUMMARY_ONLY',
        observational_support_kind='SUBSET_COVERAGE_NOT_PHYSICAL_OBSERVABILITY',
    )
    gates=admissibility(contract,has_physical_extent=False,has_observational_extent=False,
                        scale_matched=False,baseline_preregistered=False)
    files=sorted((ROOT/'data').glob('*.csv'))
    manifest={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in files}
    result={
        'schema':'dcae.reanalysis.2026-11-18.v0.2',
        'study_kind':'reanalysis_of_published_tabular_results_no_raw_observation_rediscovery',
        'source_manifest':manifest,
        'cosmos_matched_selection':comparison,
        'cosmos_boundary_support':field,
        'scsn_degradation':scsn_check,
        'dcae_cosmic_web_gate':gates,
        'can_claim_empirical_universal_growth':False,
        'can_claim_unique_physical_mechanism':False,
        'can_claim_identifiability_from_coverage_only':False,
        'new_independent_measurements':0,
    }
    path=ROOT/'results'/'REANALYSIS.json';path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    report=compute()
    print(json.dumps({
      'matched_gap':report['cosmos_matched_selection']['geometric_gap_difference_abs'],
      'matched_eps_contrasts':{k:v['walk_minus_loss_pp'] for k,v in report['cosmos_matched_selection']['tolerance_contrasts'].items()},
      'field':report['cosmos_boundary_support'],
      'gate':report['dcae_cosmic_web_gate']['state'],
      'data_sha1_verified_separately':report['scsn_degradation']['source_git_blob_sha1'],
    },indent=2,ensure_ascii=False))

import math,sys,unittest,pathlib,csv,hashlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
from dcae import scale_asymmetry, completeness_proxy, ScaleContract, admissibility
from run_analysis import compute

class DCAEAlgebra(unittest.TestCase):
  def test_reference_is_zero(self):
    x=scale_asymmetry(1200,700,1200,700)
    self.assertEqual((x.increased_0_1,x.decreased_0_1,x.direction),(0.0,0.0,'STABLE'))
  def test_identical_proportional_scaling(self):
    # An identity check, not a simulated physical observation.
    x=scale_asymmetry(2,1,20,10)
    self.assertAlmostEqual(x.delta,0,places=12)
  def test_growth_only_positive(self):
    x=scale_asymmetry(2,1,4,1)
    self.assertAlmostEqual(x.increased_0_1,0.5,places=12)
    self.assertEqual(x.decreased_0_1,0)
  def test_contraction_only_negative(self):
    x=scale_asymmetry(4,1,2,1)
    self.assertEqual(x.increased_0_1,0)
    self.assertAlmostEqual(x.decreased_0_1,0.5)
  def test_unit_invariance(self):
    x=scale_asymmetry(23.7,12.3,31,15)
    y=scale_asymmetry(23700,12300,31000,15000)
    self.assertAlmostEqual(x.delta,y.delta,places=12)
  def test_invalid_input_abstains(self):
    for v in [0,-1,float('nan'),float('inf')]:
      with self.assertRaises(ValueError): scale_asymmetry(1,1,v,1)
  def test_coverage_equivalence(self):
    for f0,f1 in [(0.7,0.9),(0.9,0.7)]:
      a=completeness_proxy(f0,f1,2)
      self.assertAlmostEqual(a.delta,math.log(f0/f1)/2)

class RealTableAudit(unittest.TestCase):
  def test_git_blob_original_scsn(self):
    # Verifies byte-for-byte transcription against actual public GitHub blob ID.
    b=(ROOT/'data'/'scsn_degradation.csv').read_bytes()
    g=hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()
    self.assertEqual(g,'54225872f65f7949527a6a8ae049868fd275e3dd')
  def test_canonical_table_intersections(self):
    r=compute()
    c=r['cosmos_matched_selection']
    self.assertEqual(c['tolerance_contrasts']['1.10']['preferred_under_this_threshold'],'independent_loss')
    self.assertEqual(c['tolerance_contrasts']['1.54']['preferred_under_this_threshold'],'walk')
    self.assertEqual(c['component_relative_rmse']['L_hd']['lower_error_mechanism'],'walk')
    self.assertEqual(c['component_relative_rmse']['f_hd']['lower_error_mechanism'],'independent_loss')
  def test_count_support_improves_on_wider_COSMOS(self):
    r=compute()['cosmos_boundary_support']
    self.assertLess(r['interior_frac_small'],r['interior_frac_large'])
    self.assertEqual(r['carrier_support_asymmetry']['increased_0_1'],0.0)
    self.assertGreater(r['carrier_support_asymmetry']['decreased_0_1'],0.0)
  def test_strict_no_promotion_to_cosmic_measurement(self):
    r=compute()
    self.assertEqual(r['dcae_cosmic_web_gate']['state'],'NOT_ESTIMABLE_DCAE')
    self.assertEqual(r['new_independent_measurements'],0)
    self.assertFalse(r['can_claim_empirical_universal_growth'])

class FKPMathematics(unittest.TestCase):
  def test_exact_volume_definition(self):
    from fkp import FkpCell
    # Algebraic unit test only; does not constitute a mock astrophysical dataset.
    c=FkpCell('a',1,8,0.1,10,0.2)
    self.assertAlmostEqual(c.information_fraction,0.25)
    self.assertAlmostEqual(c.v_eff,2.0)
  def test_input_guards(self):
    from fkp import FkpCell,cumulative_fkp
    with self.assertRaises(ValueError): FkpCell('x',1,-1,.01,10,.1)
    with self.assertRaises(ValueError): cumulative_fkp([])
    with self.assertRaises(ValueError): cumulative_fkp([FkpCell('a',1,1,.1,1,.1),FkpCell('b',2,1,.1,1,.2)])
  def test_small_nP_against_high_precision(self):
    from decimal import Decimal, localcontext
    from fkp import FkpCell
    with localcontext() as ctx:
      ctx.prec=90
      for exponent in (-5,-12,-20,-100):
        x=Decimal(10)**exponent
        reference=float((x/(1+x))**2)
        actual=FkpCell('tiny',1,1,float(x),1,.2).information_fraction
        self.assertGreater(actual,0)
        self.assertLess(abs(actual/reference-1),5e-15)
  def test_large_nP_limit(self):
    from fkp import FkpCell
    self.assertEqual(FkpCell('large',1,1,1e200,1e200,.2).information_fraction,1.0)

class InclusionProbabilityIdentity(unittest.TestCase):
  def test_union_of_two_independent_walks_by_exact_enumeration(self):
    # Exact rational verification on a three-node path, not COSMOS data.
    from fractions import Fraction as Q
    adjacency={0:(1,),1:(0,2),2:(1,)}
    walks=[((i,),Q(1,3)) for i in adjacency]
    for _ in range(2):
      walks=[(path+(j,),p/Q(len(adjacency[path[-1]])))
             for path,p in walks for j in adjacency[path[-1]]]
    for i in adjacency:
      q=sum((p for path,p in walks if i in path),Q(0))
      actual=sum((p*r for a,p in walks for b,r in walks if i in a or i in b),Q(0))
      self.assertEqual(actual,1-(1-q)**2)
    ai=sum((p for path,p in walks if 0 not in path),Q(0))
    aj=sum((p for path,p in walks if 2 not in path),Q(0))
    aij=sum((p for path,p in walks if 0 not in path and 2 not in path),Q(0))
    actual=sum((p*r for a,p in walks for b,r in walks
                if 0 in set(a+b) and 2 in set(a+b)),Q(0))
    self.assertEqual(actual,1-ai**2-aj**2+aij**2)

if __name__=='__main__': unittest.main()

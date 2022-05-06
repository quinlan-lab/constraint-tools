from snvs import fetch_SNVs
from colorize import print_string_as_info_dim

def get_M_K(window, mutations, genome, meta, number_chromosomes_min, log): 
  SNVs = fetch_SNVs(mutations, genome, window['region'], meta, number_chromosomes_min) 
  SNV_count = len(SNVs)
  singleton_count = len([SNV for SNV in SNVs if SNV['number_ALT_chromosomes'] == 1])
  if log: 
    print_string_as_info_dim(f"{singleton_count} of {SNV_count} SNVs in {window['region']} are singletons")
  return SNV_count, singleton_count

def get_M_K_trainingTime(window, mutations, genome, args): 
  return get_M_K(
    window, 
    mutations, 
    genome, 
    meta=args.__dict__, 
    number_chromosomes_min=args.number_chromosomes_min, 
    log=True)

def get_M_K_testTime(window, mutations, genome, model): 
  return get_M_K(
    window, 
    mutations, 
    genome, 
    meta=model, 
    number_chromosomes_min=model['numberChromosomesMin'], 
    log=False
  )
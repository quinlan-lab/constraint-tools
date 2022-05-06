import axios from 'axios'

let axiosConfig = {
  withCredentials: false,
}

if (process.env.VUE_APP_BASE_URL) { 
  axiosConfig['baseURL'] = process.env.VUE_APP_BASE_URL
}

console.log('axiosConfig')
console.log(axiosConfig)

const axiosInstance = axios.create(axiosConfig)

export async function getInitialPlotParameters() {
  try {
    const response = await axiosInstance.get('/api/initial-plot-parameters')
    const initialPlotParameters = response.data 
    return initialPlotParameters
  } catch (error) {
    console.error(error)
  }
}

export async function getModelParameters() {
  try {
    const response = await axiosInstance.get('/api/model-parameters')
    const modelParameters = response.data 
    return modelParameters
  } catch (error) {
    console.error(error)
  }
}

export async function getNeutralRegions(plotParameters) {
  try {
    const response = await axiosInstance.post('/api/neutral-regions', plotParameters)
    let neutralRegions = response.data['neutralRegions'] 
    console.log('neutralRegions lengths:')
    console.log(neutralRegions.map(region => region.End - region.Start))
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Functions/Arrow_functions#advanced_syntax
    neutralRegions = neutralRegions.map(region => ({
      'start': region.Start, 
      'end': region.End
    }))
    return neutralRegions
  } catch (error) {
    console.error(error)
  }
}

export async function getDistributions(payload) {
  try {
    const response = await axiosInstance.post('/api/distributions', payload)
    const distributions = response.data 
    console.log('distributions:')
    console.log(distributions)
    return distributions  
  } catch (error) {
    console.error(error)
  }
}

export async function getExpectedObservedCounts(plotParameters) {
  try {
    const response = await axiosInstance.post('/api/expected-observed-counts', plotParameters)
    const expectedObservedCounts = response.data 
    console.log('expectedObservedCounts:')
    console.log(expectedObservedCounts)
    return expectedObservedCounts  
  } catch (error) {
    console.error(error)
  }
}

function getDomain(genomeBuild) {
  if ( genomeBuild === 'hg38' ) {
    return 'https://rest.ensembl.org'
  } 
  if ( genomeBuild === 'hg37' ) {
    return 'https://grch37.rest.ensembl.org'
  } 
  throw `genome build ${genomeBuild} not supported!`
}

export async function getCanonicalTranscript(transcriptIDs, genomeBuild) {
  try { 
    // https://rest.ensembl.org/documentation/info/lookup_post
    // https://github.com/Ensembl/ensembl-rest/wiki/Getting-Started
    const response = await axiosInstance.post(`${getDomain(genomeBuild)}/lookup/id`, {
      'ids': transcriptIDs
    }, {
      params: {
        'expand': 1
      }
    })
    const transcripts = response.data 

    // https://uswest.ensembl.org/info/genome/genebuild/canonical.html
    // eslint-disable-next-line no-unused-vars
    let canonicalTranscript = Object.entries(transcripts).filter(([transcriptID, transcriptObject]) => transcriptObject.is_canonical == 1)
    if ( canonicalTranscript.length == 0 ) {
      console.log('no gene overlaps the region of interest')
      return {}
    } 
    if ( canonicalTranscript.length == 1 ) {
      // only one gene overlaps the region of interest
      canonicalTranscript = canonicalTranscript[0][1]
      return canonicalTranscript
    } 
    console.error('There is more than one canonical transcript in region.')
  } catch (error) {
    console.error(error)
  }
}

export async function getExons(region, genomeBuild) {
  try { 
    // https://rest.ensembl.org/documentation/info/overlap_region
    const response = await axios.get(`${getDomain(genomeBuild)}/overlap/region/human/${region}`, {
      params: {
        'feature': 'exon',
      }
    })
    const exons = response.data
    return exons 
  } catch (error) {
    console.error(error)
  }
}

export async function getChromosomeLength(chromosome, genomeBuild) {
  try { 
    // https://rest.ensembl.org/documentation/info/assembly_info
    const response = await axios.get(`${getDomain(genomeBuild)}/info/assembly/homo_sapiens`, {
      params: {
        'content-type': 'application/json',
      }
    })
    const contigs = response.data.top_level_region
    for ( let contig of contigs ) { 
      if ( 
        contig.coord_system === 'chromosome' && 
        contig.name === chromosome.replace('chr','') 
      ){ 
        const chromosomeLength = contig.length
        console.log(`length of chromosome ${chromosome}: ${chromosomeLength}`)
        return chromosomeLength
      }
    }
    throw `length of chromosome ${chromosome} cannot be obtained from Ensemble API!`
  } catch (error) {
    console.error(error)
  }
}
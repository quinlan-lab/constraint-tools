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

export async function getCanonicalTranscript(transcriptIDs) {
  try { 
    // https://rest.ensembl.org/documentation/info/lookup_post
    // https://github.com/Ensembl/ensembl-rest/wiki/Getting-Started
    const response = await axiosInstance.post('https://rest.ensembl.org/lookup/id', {
      'ids': transcriptIDs
    }, {
      params: {
        'expand': 1
      }
    })
    const transcripts = response.data 

    // http://uswest.ensembl.org/info/genome/genebuild/canonical.html
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

export async function getExons(region) {
  try { 
    // https://rest.ensembl.org/documentation/info/overlap_region
    const response = await axios.get(`http://rest.ensembl.org/overlap/region/human/${region}`, {
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
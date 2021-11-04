import axios from 'axios'

const axiosInstance = axios.create({
  // baseURL: 'http://localhost:5000', // toggle this; required for development of vue app
  withCredentials: false,
})

export async function getInitialPlotParameters() {
  try {
    const response = await axiosInstance.get('/api/initial-plot-parameters')
    const initialPlotParameters = response.data 
    return initialPlotParameters
  } catch (error) {
    console.error(error)
  }
}

export async function getMutationCounts(plotParameters) {
  try {
    const response = await axiosInstance.post('/api/mutation-counts', plotParameters)
    const mutationCounts = response.data 
    console.log('mutationCounts:')
    console.log(mutationCounts)
    return mutationCounts  
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
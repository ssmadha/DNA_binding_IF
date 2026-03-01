const [
    asif,
    asif_melted,
    gene_cluster,
    tissue_cluster,
    ppi_edges_byGene,
    ppi_edges_string_byGene,
    hi_union
] = await Promise.all([
    d3.csv("./asif/ASIF.csv"),
    d3.csv("./asif/ASIF_melted.csv"),
    d3.csv("./asif/gene_clusters.csv"),
    d3.csv("./asif/tissue_clusters.csv"),
    d3.csv("./asif/ppi_edges_byGene.csv"),
    d3.csv("./asif/ppi_edges_string_byGene.csv"),
    d3.csv("./asif/HI-union_TFonly.csv")
]);

//import * as d3 from 'd3';

//read data so that every column but the first is numeric
asif.forEach((d) => {
  Object.keys(d).forEach((key) => {
    if (key !== 'Gene ID' && key !== 'Gene Name') {
      d[key] = +d[key];
    }
  });
});

for (const d of asif_melted) {
  d.ASIF = +d.ASIF;
  d.Expression = +d.Expression;
}

//asif_melted.sort((a, b) => {
//  return (
//    Object.keys(gene_cluster[0]).indexOf(b['Gene ID']) -
//      Object.keys(gene_cluster[0]).indexOf(a['Gene ID']) ||
//    Object.keys(tissue_cluster[0]).indexOf(a['Tissue']) -
//      Object.keys(tissue_cluster[0]).indexOf(b['Tissue'])
//  );
//});

// below is optimized version of above from ChatGPT
//
// Preprocess keys
const geneKeys = Object.keys(gene_cluster[0]);
const tissueKeys = Object.keys(tissue_cluster[0]);

// Create maps for faster lookup
const geneKeyIndexes = {};
geneKeys.forEach(
  (key, index) => (geneKeyIndexes[key] = index),
);

const tissueKeyIndexes = {};
tissueKeys.forEach(
  (key, index) => (tissueKeyIndexes[key] = index),
);

// Sort the array
asif_melted.sort((a, b) => {
  return (
    geneKeyIndexes[b['Gene ID']] -
      geneKeyIndexes[a['Gene ID']] ||
    tissueKeyIndexes[a['Tissue']] -
      tissueKeyIndexes[b['Tissue']]
  );
});

for (const key of Object.keys(gene_cluster[0])) {
  gene_cluster[0][key] = +gene_cluster[0][key];
}

for (const key of Object.keys(tissue_cluster[0])) {
  tissue_cluster[0][key] = +tissue_cluster[0][key];
}

for (const d of ppi_edges_byGene) {
  d.coverage = +d.coverage;
}

for (const d of ppi_edges_string_byGene) {
  d.coverage = +d.coverage;
}

for (const d of hi_union) {
  d.coverage = +d.coverage;
}
//import data from this URL: https://raw.githubusercontent.com/ssmadha/DNA_binding_IF/main/ASIF_transcript_melted_no0.csv
//import axios from 'axios';

const url =
  'https://raw.githubusercontent.com/ssmadha/DNA_binding_IF/main/ASIF_transcript_melted_no0.csv';

let transcript_data;
async function fetchData() {
  try {
    // Use await to wait for the promise to resolve
    const data = await d3.csv(url);
    // Process the data or perform any other actions
    for (const d of data) {
      d.ASIF = +d.ASIF;
      d.Expression = +d.Expression;
      d['Percent Expression'] = +d['Percent Expression'];
    }
    return data; // Return the parsed CSV data
  } catch (error) {
    // Handle any errors that occur during fetching or parsing
    console.error('Error loading the CSV file:', error);
  }
}

// Call the fetchData function to fetch and parse CSV data
fetchData().then((data) => {
  // Store the data in a variable or perform further actions
  //console.log(data); // Output the parsed CSV data to the console
  transcript_data = data;
  // You can assign the data to a variable here if needed
});
export const main = (container) => {
  const fontSize = 20;

  //const json = JSON.stringify(transcript_data[0], null, 2);

  // container.innerHTML = `
  //   <pre style="font-size: ${fontSize}px;">${json}</pre>
  // `;
};

export {
  asif,
  asif_melted,
  gene_cluster,
  tissue_cluster,
  ppi_edges_byGene,
  ppi_edges_string_byGene,
  hi_union,
  transcript_data,
};

//import * as d3 from 'd3';

export const infoView = (
  selection,
  { nodes, data, subData, tissue },
) => {
  //console.log(nodes);

  //find the row in data whose ID column matches node.id
  const matchedRows = nodes
    ? data
        .filter((d) => {
          return (
            nodes.some((node) => {
              return d['Gene ID'] === node.id;
            }) && d.Tissue === tissue
          );
        })
        .map((d) => {
          return {
            Name: d['Gene Name'],
            ID: d['Gene ID'],
            ASIF: d.ASIF.toFixed(3),
            Expression: d.Expression.toFixed(2),
          };
        })
    : [];
  const table = selection
    .selectAll('table.info-view')
    .data([null])
    .join('table')
    .attr('class', 'table info-view');
  const geneHeader = ['Name', 'ID', 'ASIF', 'Expression'];

  // Check if the thead element exists
  if (table.select('thead').empty()) {
    // If not, add the header to the table
    table
      .append('thead')
      .append('tr')
      .selectAll('th')
      .data(geneHeader)
      .join('th')
      .text((d) => d);
  }

  // Add row based on matchedRow
  const tbody = table.select('tbody').empty()
    ? table
        .append('tbody')
        .attr('class', 'table-group-divider')
    : table.select('tbody');
  if (!subData) {
    return;
  }
  const subRowsData = subData
    .filter(
      (d) =>
        d.Tissue === tissue &&
        d['Gene'] === matchedRows[0].ID,
    )
    .map((d) => {
      return {
        ID: d['Transcript'],
        ASIF: d.ASIF.toFixed(3),
        Expression: d.Expression.toFixed(2),
        'Percent Expression':
          (d['Percent Expression'] * 100).toFixed(2) + '%',
      };
    })
    .sort((a, b) => +a.Expression < +b.Expression);

  tbody.selectAll('tr').remove();

  // Join the data with the table rows
  tbody
    .selectAll('tr.gene')
    .data(matchedRows, (d) => d.ID)
    .join(
      (enter) => {
        const geneRow = enter
          .append('tr')
          .attr('class', 'gene');

        geneRow
          .selectAll('td')
          .data((d) => Object.values(d))
          .join('td')
          .text((d) => d);

        // geneRow
        //   .append('td')
        //   .append('button')
        //   .text('Show Transcripts');

        const subTable = enter
          .append('tr')
          .attr('class', 'sub-gene')
          .append('td')
          .attr('colspan', 4)
          .append('table')
          .attr('class', 'transcript table');

        subTable
          .append('thead')
          .append('tr')
          .selectAll('th')
          .data(Object.keys(subRowsData[0]))
          .join('th')
          .text((d) => d);

        subTable
          .append('tbody')
          .attr('class', 'table-divider')
          .selectAll('tr')
          .data(subRowsData, (d) => d.ID)
          .join('tr')
          .selectAll('td')
          .data((d) => Object.values(d))
          .join('td')
          .text((d) => d);
      },
      (update) => {
        update;
      },
      (exit) => {
        exit.remove();
      },
    );
};

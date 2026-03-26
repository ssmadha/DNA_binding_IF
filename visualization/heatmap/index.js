import { select } from 'd3';
import { viz } from './viz';
import {
  asif_melted,
  gene_cluster,
  tissue_cluster,
} from '@ssmadha/asif';
import { observeResize } from '@ssmadha/responsive-axes';

const tissueColumn = 'Tissue';
const nameColumn = 'Gene Name';
const idColumn = 'Gene ID';
const valueColumn = 'ASIF';

export const main = (container, { state, setState }) => {
  const dimensions = observeResize({
    state,
    setState,
    container,
  });

  if (dimensions === null) return;

  const { width, height } = dimensions;

  const svg = select(container)
    .selectAll('svg')
    .data([null])
    .join('svg')
    .attr('width', width)
    .attr('height', height);

  const vizState = {
    ...state,
    data: asif_melted,
    xLabels: (d) => d[tissueColumn],
    xRotation: 45,
    xAxisLabelText: 'Tissue',
    xAxisLabelOffset: 80,
    yLabels: (d) => d[nameColumn],
    yTickSize: 0,
    yAxisLabelText: '',
    yAxisLabelOffset: 17,
    heatValue: (d) => d[valueColumn],
    innerRectFill: '#E8E8E8',
    padding: 0 / 100,
    marginTop: 20,
    marginBottom: 95,
    marginLeft: 50,
    marginRight: 100,
    width,
    height,
    xAnnotation: tissue_cluster[0],
    xAnnoPadding: 4,
    xAnnoHeight: 10,
    yAnnotation: gene_cluster[0],
    yAnnoLabels: (d) => d[idColumn],
    yAnnoPadding: 4,
    yAnnoWidth: 10,
    legendLabel: 'ASIF',
    legendX: width - 70,
    legendY: 30,
    legendPadding: 20,
    legendHeight: 250,
    legendWidth: 20,
    legendTickNumber: 8,
    legendColorNumber: 100,
  };

  viz(svg, {
    state: vizState,
    setState,
  });
};

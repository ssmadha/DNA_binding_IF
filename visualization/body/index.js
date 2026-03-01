//import { select } from 'd3';
import { body } from './body.js';
import { bodyData } from './bodyData.js';
import { organDict } from './organDict.js';

export const main = (container, { state, setState }) => {
  const width = container.clientWidth;
  const height = container.clientHeight;

  const { clickedOrgan1, clickedOrgan2, hoveredOrgan } =
    state;

  const setClickedOrgan = (
    clickedOrgan1,
    clickedOrgan2,
  ) => {
    setState((state) => ({
      ...state,
      clickedOrgan1,
      clickedOrgan2,
    }));
  };

  const setHoveredOrgan = (hoveredOrgan) => {
    setState((state) => ({
      ...state,
      hoveredOrgan,
    }));
  };

  const svg = select(container)
    .selectAll('svg.body')
    .data([null])
    .join('svg')
    .attr('class', 'body')
    .attr('viewbox', '0 0 1000 1000')
    .attr('width', '100%')
    .attr('height', '100%')
    //.attr('preserveAspectRatio', 'none')
    .style('background', '#FFFFFF');

  svg.call(body, {
    bodyData,
    organDict,
    legendLabel: 'Organ',
    legendLabelFontSize: 22,
    legendX: 30,
    legendY: 36,
    legendPadding: 17,
    legendFontSize: 12,
    tickSpacing: 11,
    tickPadding: 0,
    clickedOrgan1,
    clickedOrgan2,
    setClickedOrgan,
    hoveredOrgan,
    setHoveredOrgan,
  });
};
export { bodyData, organDict };

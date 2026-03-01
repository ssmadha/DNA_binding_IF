//import * as d3 from 'd3';

export const simulate = ({
  nodes,
  links,
  linkIDs = (d) => d.id,
  netGraphWidth,
  netGraphHeight,
  simMinX,
  simMaxX,
  simMinY,
  simMaxY,
  linkDistance,
  linkStrength,
  linkIterations,
  manyBodyStrength,
  manyBodyMinDist,
  manyBodyMaxDist,
  collideRadius,
  centerX,
  centerY,
  centerStrength,
  strengthXY,
  velocityDecay,
}) => {
  const sim = d3
    .forceSimulation(nodes)
    .force(
      'link',
      d3
        .forceLink(links)
        .id(linkIDs)
        .distance(linkDistance)
        .strength(linkStrength)
        .iterations(linkIterations),
    )
    .force(
      'charge',
      d3
        .forceManyBody()
        .strength(manyBodyStrength)
        .distanceMin(manyBodyMinDist)
        .distanceMax(manyBodyMaxDist),
    )
    .force(
      'collide',
      d3.forceCollide().radius(collideRadius),
    )
    // .force(
    //   'center',
    //   d3
    //     .forceCenter(centerX, centerY)
    //     .strength(centerStrength),
    // )
    .force(
      'x',
      d3
        .forceX(centerX)
        .strength((d) =>
          d.x < simMinX && d.x > simMaxX
            ? strengthXY
            : strengthXY / 2,
        ),
    )
    .force(
      'y',
      d3
        .forceY(centerY)
        .strength((d) =>
          d.y < simMinY && d.y > simMaxY
            ? strengthXY
            : strengthXY / 2,
        ),
    )
    .velocityDecay(velocityDecay)
    .stop();
  // .on('tick', () =>
  //   graph(graphGroup, {
  //     nodes,
  //     links,
  //     transform,
  //     setState,
  //   }),
  //);

  sim.tick(
    Math.ceil(
      Math.log(sim.alphaMin()) /
        Math.log(1 - sim.alphaDecay()),
    ),
  );
};

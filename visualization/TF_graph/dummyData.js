export const data = [
  { 'Gene Name': 'A', Expression: 20, asif: 0.1 },
  { 'Gene Name': 'B', Expression: 13, asif: 0.9 },
  { 'Gene Name': 'C', Expression: 20, asif: 0.45 },
  { 'Gene Name': 'D', Expression: 12, asif: 0.1 },
  { 'Gene Name': 'E', Expression: 12, asif: 0.95 },
  { 'Gene Name': 'F', Expression: 19, asif: 0.75 },
  { 'Gene Name': 'G', Expression: 16, asif: 0.25 },
  { 'Gene Name': 'H', Expression: 17, asif: 0.0 },
];

export const edges = [
   { source: 'A', target: 'B', coverage: 0.0 },
   { source: 'C', target: 'D', coverage: 1.0 },
   { source: 'E', target: 'G', coverage: 0.5 },
   { source: 'G', target: 'E', coverage: 0.75 },
   { source: 'G', target: 'A', coverage: 0.25 },
   { source: 'F', target: 'A', coverage: 0.0 },
   { source: 'H', target: 'H', coverage: 1.0 },
];

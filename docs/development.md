# Development

## Style conventions

It would be handy to be able to differentiate between modules, subworkflow, and workflows at a glance.
To make that happen, use the following naming conventions that imply a hierarchy:

- Modules: lower_snake_case
  - Example: `samtools_sort`
- Subworkflow: Pascal_Snake_Case
  - Example: `Prepare_Refs`
- Workflow: UPPER_SNAKE_CASE
  - Example: `PROCESS_READS`

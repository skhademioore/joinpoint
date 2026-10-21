# Security Policy

## Supported Versions

Security fixes are considered for the active development branch and the
most recent release.

| Version   | Supported |
|-----------|-----------|
| `0.1.x`   | Yes       |
| `< 0.1.0` | No        |

## Reporting a Vulnerability

Please do not report security vulnerabilities in public GitHub issues.

Report suspected vulnerabilities by email to
`sahar.khademi975@gmail.com`, or use GitHub private vulnerability
reporting if it is enabled for the repository.

Include as much detail as you can:

- affected version or commit,
- operating system and R version,
- minimal reproduction steps,
- whether the issue affects local execution, package installation, the
  compiled C++ code, file handling, or documentation output,
- any proof-of-concept code.

## Response Process

The maintainer will acknowledge reports as soon as practical,
investigate the issue, and coordinate a fix if the report is confirmed.
Public disclosure should wait until a fix or mitigation is available,
unless there is an overriding user safety reason to disclose earlier.

## Scope

Security-sensitive areas include:

- arbitrary file writes or path traversal (for example through
  [`jp_export()`](https://skhademioore.github.io/joinpoint/reference/jp_export.md)),
- unsafe execution during package installation or compilation,
- exposure of private data in examples, tests, or logs,
- documentation content that could mislead users into unsafe use.

Methodological concerns that do not create a software-security
vulnerability should be opened as normal GitHub issues instead.

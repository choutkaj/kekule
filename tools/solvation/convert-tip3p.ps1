# Offline conversion of the pinned OpenMM TIP3P PDB. No simulation or fitting.
param([Parameter(Mandatory=$true)][string]$SourcePdb)
$ErrorActionPreference = 'Stop'
$expected = '14FD37900D627C0E258D6086A14C6084E4BEC1422E9D20FCCFC83C3F814FD7EE'
if ((Get-FileHash -LiteralPath $SourcePdb -Algorithm SHA256).Hash -ne $expected) {
    throw 'TIP3P source checksum mismatch'
}
$outputPath = Join-Path $PSScriptRoot '../../crates/kekule/src/structure/solvation/data/tip3p.bin'
$bytes = [IO.MemoryStream]::new()
$writer = [IO.BinaryWriter]::new($bytes)
$count = 0
foreach ($line in [IO.File]::ReadLines((Resolve-Path -LiteralPath $SourcePdb))) {
    if (!$line.StartsWith('ATOM  ')) { continue }
    $expectedName = @('O', 'H1', 'H2')[$count % 3]
    if ($line.Substring(12,4).Trim() -ne $expectedName) { throw 'Unexpected water atom ordering' }
    foreach ($offset in @(30,38,46)) {
        $value = [decimal]::Parse($line.Substring($offset,8), [Globalization.CultureInfo]::InvariantCulture)
        $writer.Write([int32]($value * 1000))
    }
    $count++
}
if ($count -ne 2685) { throw 'Unexpected atom count' }
$writer.Flush()
[IO.File]::WriteAllBytes($outputPath, $bytes.ToArray())
$writer.Dispose()
$bytes.Dispose()
(Get-FileHash -LiteralPath $outputPath -Algorithm SHA256).Hash

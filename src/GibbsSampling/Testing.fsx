

#r @"..\GibbsSampling\bin\Release\netstandard2.0\FSharp.Plotly.dll"
#r @"..\GibbsSampling\bin\Release\netstandard2.0\FSharpAux.dll"
#r @"..\GibbsSampling\bin\Release\netstandard2.0\FSharpAux.IO.dll"
#r @"..\GibbsSampling\bin\Release\netstandard2.0\BioFSharp.dll"
#r @"..\GibbsSampling\bin\Release\netstandard2.0\BioFSharp.Stats.dll"
#r @"..\GibbsSampling\bin\Release\netstandard2.0\BioFSharp.IO.dll"
#r @"..\GibbsSampling\bin\Release\netstandard2.0\GibbsSampling.dll"


open BioFSharp
open FSharpAux
open GibbsSampling
open SiteSampler

#time

let tests =
    [|
        "GTGGCTGCACCACGTGTATGC"
        "ACATCGCATCACGTGACCAGT"
        "CCTCGCACGTGGTGGTACAGT"
        "CTCGTTAGGACCATCACGTGA"
    |]

let profileSequenceForTests =
    [
        "CACGTG"
        "CACGTG"
        "CACGTG"
        "CACGTG"
    ]

let bioTests =
    tests
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let bioTestsWithMultipleSamples = 
    [|
        "GTGGCTGCACCACGTGTATGCCACGTG"
        "ACATCGCATCACGTGACCAGTTAGTTG"
        "CCTCGCACGTGGTGGTACAGTCGTACG"
        "GCATAAAGGACCATCACGTGAAGCTGC"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTT"
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let bioTestsII =
    [|
        "GTAAGTACAGAAAGCCACAGAGTACCATCTAGGAAATTAACATTATACTAACTTTCTACATCGTTGATACTTATGCGTATACATTCATATA"
        "AGACAGAGTCTAAAGATTGCATTACAAGAAAAAAGTTCTCATTACTAACAAGCAAAATGTTTTGTTTCTCCTTTTA"
        "GTATGTTCATGTCTCATTCTCCTTTTCGGCTCCGTTTAGGTGATAAACGTACTATATTGTGAAAGATTATTTACTAACGACACATTGAAG*"
        "GCATGTGTGCTGCCCAAGTTGAGAAGAGATACTAACAAAATGACCGCGGCTCTCAAAAATAATTGACGAGCTTACGGTGATACGCTTACCG"
        "GTATGTTTGACGAGAATTGCTAGTGTGCGGGAAACTTTGCTACCTTTTTTGGTGCGATGCAACAGGTTACTAATATGTAATACTTCAG"
        "TTTCAAGATTAACCACATCTGCTAACTTTCTCCCTATGCTTTTACTAACAAAATTATTCTCACTCCCCGATATTGA"
        "GTAAGTATCCAGATTTTACTTCATATATTTGCCTTTTTCTGTGCTCCGACTTACTAACATTGTATTCTCCCCTTCTTCATTTTAG"
        "GTATGCATAGGCAATAACTTCGGCCTCATACTCAAAGAACACGTTTACTAACATAACTTATTTACATAG"
        "GTATGTAGTAGGGAAATATATCAAAGGAACAAAATGAAAGCTATGTGATTCCGTAATTTACGAAGGCAAATTACTAACATTGAAATACGGG"
        "GTATGTTACTATTTGGAGTTTCATGAGGCTTTTCCCGCCGTAGATCGAACCCAATCTTACTAACAGAGAAAGGGCTTTTTCCCGACCATCA"
        "TATGTAATGATATATTATGAAGTAAGTTCCCCAAAGCCAATTAACTAACCGAATTTTAATCTGCACTCATCATTAG"
        "GTATGTTCATAATGATTTACATCGGAATTCCCTTTGATACAAGAAAACTAACGGGTATCGTACATCAATTTTTGAAAAAAGTCAAGTACTA"
        "GTATGTATATTTTTGACTTTTTGAGTCTCAACTACCGAAGAGAAATAAACTACTAACGTACTTTAATATTTATAG"
        "TTTCGACGCGAATAGACTTTTTCCTTCTTACAGAACGATAATAACTAACATGACTTTAACAG"
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let profileSequenceForTestsII =
    [|"TACTAAC"; "TACTAAT"; "AACTAAC"|]

/// Contains a Set of all DNA bases.
let dnaBases = 
    [|Nucleotides.Nucleotide.A; Nucleotides.Nucleotide.T; Nucleotides.Nucleotide.G; Nucleotides.Nucleotide.C; Nucleotides.Nucleotide.Gap|]

/// Contains a Set of all amino acids.
let aminoAcids  =   
    [|
        AminoAcidSymbols.AminoAcidSymbol.Ala; AminoAcidSymbols.AminoAcidSymbol.Arg; AminoAcidSymbols.AminoAcidSymbol.Asn; 
        AminoAcidSymbols.AminoAcidSymbol.Asp; AminoAcidSymbols.AminoAcidSymbol.Asx; AminoAcidSymbols.AminoAcidSymbol.Cys;
        AminoAcidSymbols.AminoAcidSymbol.Xle; AminoAcidSymbols.AminoAcidSymbol.Gln; AminoAcidSymbols.AminoAcidSymbol.Glu;
        AminoAcidSymbols.AminoAcidSymbol.Glx; AminoAcidSymbols.AminoAcidSymbol.Gly; AminoAcidSymbols.AminoAcidSymbol.His;
        AminoAcidSymbols.AminoAcidSymbol.Ile; AminoAcidSymbols.AminoAcidSymbol.Leu; AminoAcidSymbols.AminoAcidSymbol.Lys;
        AminoAcidSymbols.AminoAcidSymbol.Met; AminoAcidSymbols.AminoAcidSymbol.Phe; AminoAcidSymbols.AminoAcidSymbol.Pro;
        AminoAcidSymbols.AminoAcidSymbol.Pyl; AminoAcidSymbols.AminoAcidSymbol.Sel; AminoAcidSymbols.AminoAcidSymbol.Ser;
        AminoAcidSymbols.AminoAcidSymbol.Thr; AminoAcidSymbols.AminoAcidSymbol.Trp; AminoAcidSymbols.AminoAcidSymbol.Val;
    |]

let testI = Array.init 10000 (fun _ -> getmotifsWithBestInformationContent 10000 6 0.0001 dnaBases bioTests)

testI
|> Array.countBy (fun items -> items |> Array.map (fun item -> snd item))
|> Array.sortByDescending (fun (_, i) -> i)

let testII = Array.init 1 (fun _ -> getmotifsWithBestInformationContent 1 7 0.0001 dnaBases bioTestsII)

testII
|> Array.countBy (fun items -> items |> Array.map (fun item -> snd item))
|> Array.sortByDescending (fun (_, i) -> i)


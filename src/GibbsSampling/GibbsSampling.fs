namespace GibbsSampling


open System
open FSharpAux
open FSharpAux.IO
open BioFSharp


///Consists of functions to work with CompositeVectors.
module CompositeVector =

    /// One dimensional array with fixed positions for each element.
    type CompositeVector<'a, 'value when 'a :> IBioItem>internal (array:'value []) =
    
        let getIndex (key:'a) =
            (BioItem.symbol key |> int) - 42

        new () =
            let arr:'value [] = Array.zeroCreate 49
            new CompositeVector<_, 'value>(arr)

        member this.Array = 
            if (obj.ReferenceEquals(array, null)) then
                raise (ArgumentNullException("array"))
            array
            
        member this.Item
            with get i       = array.[getIndex i]
            and  set i value = array.[getIndex i] <- value

    /// One dimensional array with fixed positions for each element.
    /// Use to track frequency of elements independent of position in source.
    type FrequencyCompositeVector internal (array:int []) =

        inherit CompositeVector<IBioItem, int>(array)

        new() = 
            let arr:'value [] = Array.zeroCreate 49
            new FrequencyCompositeVector(arr)

    /// Increase counter of position by 1.
    let increaseInPlaceFCV (bioItem:#IBioItem) (frequencyCompositeVector:FrequencyCompositeVector) = 
        frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + 1
        frequencyCompositeVector

    /// Increase counter of position by n.
    let increaseInPlaceFCVBy (bioItem:'a) (frequencyCompositeVector:FrequencyCompositeVector) n = 
        frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + n
        frequencyCompositeVector

    ///// Create a FrequencyCompositeVector based on BioArrays and exclude the specified segments.
    //let createFCVBy (motiveLength:int) (resSources:BioArray.BioArray<#IBioItem>[]) (motivePositions:int[]) =
    //    let backGroundCounts = new FrequencyCompositeVector()    
    //    Array.map2 (fun (resSource:BioArray.BioArray<#IBioItem>) (position:int) -> Array.append resSource.[0..(position-1)] resSource.[(position+motiveLength)..]) resSources motivePositions
    //    |> Array.concat
    //    |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

    /// Create a FrequencyCompositeVector based on BioArrays and exclude the specified segments.
    let createFCVOf (resSources:BioArray.BioArray<#IBioItem>) =
        let backGroundCounts = new FrequencyCompositeVector()   
        Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts resSources

    /// Create new FrequencyCompositeVector based on the sum of the positions of an array of FrequencyVectors.
    let fuseFrequencyVectors (alphabet:#IBioItem[]) (bfVectors:FrequencyCompositeVector[]) =
        let backgroundFrequencyVector = new FrequencyCompositeVector()
        for fcVector in bfVectors do
            for bioItem in alphabet do
                backgroundFrequencyVector.[bioItem] <- backgroundFrequencyVector.[bioItem] + (fcVector.[bioItem])
        backgroundFrequencyVector

    /// Create FrequencyCompositeVector based on BioArray and excludes the specified segment.
    let createFCVWithout (motiveLength:int) (position:int) (resSource:BioArray.BioArray<#IBioItem>) =
        let backGroundCounts = new FrequencyCompositeVector()
        Array.append resSource.[0..(position-1)] resSource.[(position+motiveLength)..]
        |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

    /// Create FrequencyCompositeVector based on BioArray.
    let increaseInPlaceFCVOf (resSources:BioArray.BioArray<#IBioItem>) (backGroundCounts:FrequencyCompositeVector) =   
        resSources
        |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

    /// Subtracts the amount of elements in the given source from the FrequencyCompositeVector.
    let substractSegmentCountsFrom (source:BioArray.BioArray<#IBioItem>) (fcVector:FrequencyCompositeVector) =
        let bfVec = new FrequencyCompositeVector(fcVector.Array)
        for bioItem in source do
            bfVec.[bioItem] <- (if fcVector.[bioItem] - 1 > 0 then fcVector.[bioItem] - 1 else 0)
        bfVec
    
    type ProbabilityCompositeVector internal (array:float []) =
        
        inherit CompositeVector<IBioItem,float>(array)

        new() = 
            let arr:'value [] = Array.zeroCreate 49
            new ProbabilityCompositeVector(arr)

    /// Increase counter of position by 1.
    let increaseInPlacePCV (bioItem:'a) (backGroundProbabilityVector:ProbabilityCompositeVector) = 
        backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + 1.
        backGroundProbabilityVector

    /// Increase counter of position by n.
    let increaseInPlacePCVBy (bioItem:'a) (backGroundProbabilityVector:ProbabilityCompositeVector) n = 
        backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + n
        backGroundProbabilityVector

    /// Create a ProbabilityCompositeVector based on a FrequencyCompositeVector by replacing the integers with floats.
    let createPCVOf (caArray:FrequencyCompositeVector) =
        caArray.Array
        |> Array.map (fun item -> float item)
        |> fun item -> new ProbabilityCompositeVector(item)

    /// Create normalized ProbabilityCompositeVector based on FrequencyCompositeVector.
    let createNormalizedPCVOfFCV (alphabet:#IBioItem[]) (pseudoCount:float) (frequencyCompositeVector:FrequencyCompositeVector) =
        let backGroundProbabilityVector = createPCVOf frequencyCompositeVector
        let sum = (float (Array.sum frequencyCompositeVector.Array)) + ((float alphabet.Length) * pseudoCount)
        for item in alphabet do
            backGroundProbabilityVector.[item] <- (backGroundProbabilityVector.[item] + pseudoCount)/sum
        backGroundProbabilityVector

    /// Calculate the score of the given sequence based on the probabilities of the ProbabilityCompositeVector.
    let calculateSegmentScoreBy (pcv:ProbabilityCompositeVector) (bioItems:BioArray.BioArray<#IBioItem>) =
        Array.fold (fun (value:float) (bios:#IBioItem) -> value * (pcv.[bios])) 1. bioItems

module PositionMatrix =

    ///Checks whether all elements in the list have a wider distance than width or not.
    let internal ceckForDistance (width:int) (items:int list) =
        if items.Length <= 0 || items.Length = 1 then true
        else
            let rec loop n i =
                if n = items.Length-1 then true
                else
                    if i >= items.Length then loop (n+1) (n+2)
                    else
                        if Operators.abs(items.[n]-items.[i]) > width then
                            loop n (i+1)
                        else false
            loop 0 1

    /// Get an integer which is between 0 and the length of the sequence - segmentLength
    let internal getRandomNumberInSequence (segmentLength:int) (source:'a[]) =
        let rnd = System.Random()
        Array.init 1 (fun _ -> rnd.Next(0, source.Length-segmentLength+1))
        |> Array.head

    /// Create a specific sub sequence of the source sequence based on the given length and starting position. Do not forget, sequences start counting at 0!
    let internal getSegment (subsequenceLength:int) (source:'a[]) (startPoint:int) =
        source 
        |> Array.skip startPoint 
        |> Array.take subsequenceLength 
        |> Array.ofSeq

    ///Finds the best information content in an array of arrays of PWMSs and positions.
    let getBestInformationContent (item:((float*int)[])[]) =
        let rec loop (n:int) (bestPWMS:(float*int)[]) =
            if n = item.Length then bestPWMS
            else
                let informationContentItem =
                    item.[n]
                    |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                let informationContentBestPWMS =
                    bestPWMS
                    |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                if informationContentItem > informationContentBestPWMS then
                    loop (n + 1) item.[n]
                else
                    loop (n + 1) bestPWMS
        loop 0 [||]

    /// Matrix with fixed positions for nucleotides and amino acids.
    type BaseMatrix<'a, 'value when 'a :> IBioItem>internal(matrix:'value [,]) =

        let getRowArray2DIndex (key:'a) =
            (BioItem.symbol key |> int) - 42

        new (rowLength:int) =
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new BaseMatrix<_, 'value>(arr)

        member this.Matrix = 
            if (obj.ReferenceEquals(matrix, null)) then
                raise (ArgumentNullException("array"))
            matrix
            
        member this.Item
            with get (column, row)       = matrix.[getRowArray2DIndex column, row]
            and  set (column, row) value = matrix.[getRowArray2DIndex column, row] <- value

    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0.
    type PositionFrequencyMatrix internal(matrix:int [,]) =

        inherit BaseMatrix<IBioItem, int>(matrix)

        new (rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionFrequencyMatrix(arr)

    /// Increase counter of PositionFrequencyMatrix at fixed position by 1.
    let increaseInPlacePFM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionFrequencyMatrix:PositionFrequencyMatrix) = 
        positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + 1
        positionFrequencyMatrix

    /// Increase counter of PositionFrequencyMatrix at fixed position by n.
    let increaseInPlacePFMBy (pos:int) (bioItem:'a when 'a :> IBioItem) (positionFrequencyMatrix:PositionFrequencyMatrix) n = 
        positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + n
        positionFrequencyMatrix

    /// Create PositionFrequencyMatrix based on BioArray.
    let createPFMOf (source:BioArray.BioArray<#IBioItem>) =
        let positionFrequencyMatrix = new PositionFrequencyMatrix(source.Length)
        source
        |> Array.fold (fun (row, cm) column -> row + 1, increaseInPlacePFM row column cm) (0, positionFrequencyMatrix) |> ignore
        positionFrequencyMatrix

    /// Create new PositionFrequencyMatrix based on the sum of the positions of an array of Countmatrices.
    let fusePositionFrequencyMatrices (motiveLength:int) (countMatrices:PositionFrequencyMatrix[]) =
        let positionFrequencyMatrix = new PositionFrequencyMatrix(motiveLength)
        if countMatrices.Length > 0 then
            for cMatrix in countMatrices do
                for column = 0 to (Array2D.length1 cMatrix.Matrix)-1 do
                        for row = 0 to (Array2D.length2 cMatrix.Matrix)-1 do
                            positionFrequencyMatrix.Matrix.[column, row] <- positionFrequencyMatrix.Matrix.[column, row] + (cMatrix.Matrix.[column, row])
            positionFrequencyMatrix
        else positionFrequencyMatrix 
    
    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0. .
    type PositionProbabilityMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionProbabilityMatrix(arr)

    /// Increase counter of PositionProbabilityMatrix at fixed position by 1.
    let increaseInPlacePPM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionProbabilityMatrix:PositionProbabilityMatrix) = 
        positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + 1.
        positionProbabilityMatrix

    /// Increase counter of PositionProbabilityMatrix at fixed position by n.
    let increaseInPlacePPMBy (pos:int) (bioItem:'a when 'a :> IBioItem) (positionProbabilityMatrix:PositionProbabilityMatrix) n = 
        positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + n
        positionProbabilityMatrix

    /// Create new PositionWeightMatrix based on existing PositionFrequencyMatrix. 
    /// The counts of each position of each element are transformed to floats.
    let createPPMOf (positionFrequencyMatrix:PositionFrequencyMatrix) =
        positionFrequencyMatrix.Matrix |> Array2D.map (fun item -> float item)
        |> fun item -> new PositionProbabilityMatrix(item)

    // Create new PositionWeightMatrix based on existing PositionFrequencyMatrix. 
    /// The counts of each position of each element are transformed to floats.
    let normalizePPM (sourceCount:int) (alphabet:#IBioItem[]) (pseudoCount:float) (positionFrequencyMatrix:PositionProbabilityMatrix) =
        let positionProbabilityMatrix = new PositionProbabilityMatrix(positionFrequencyMatrix.Matrix)
        let sum = (float sourceCount) + ((float alphabet.Length) * pseudoCount)
        for item in alphabet do
            for position = 0 to (Array2D.length2 positionProbabilityMatrix.Matrix) - 1 do
            positionProbabilityMatrix.[item, position] <- (positionProbabilityMatrix.[item, position] + pseudoCount)/sum
        positionProbabilityMatrix
    
    type PositionWeightMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionWeightMatrix(arr)

    /// Increase counter of PositionWeightMatrix at fixed position by 1.
    let increaseInPlacePWM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionWeightMatrix:PositionWeightMatrix) = 
        positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + 1.
        positionWeightMatrix

    // Increase counter of PositionWeightMatrix at fixed position by n.
    let increaseInPlacePWMBy (pos:int) (bioItem:'a when 'a :> IBioItem) (positionWeightMatrix:PositionWeightMatrix) n = 
        positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + n
        positionWeightMatrix

    /// Create PositionWeightMatrix based on ProbabilityCompositeVector and PositionProbabilityMatrix.
    let createPositionWeightMatrix (alphabet:#IBioItem[]) (pcv:CompositeVector.ProbabilityCompositeVector) (ppMatrix:PositionProbabilityMatrix) =
        let pwMatrix = new PositionWeightMatrix(Array2D.length2 ppMatrix.Matrix)
        for item in alphabet do
            for position=0 to (Array2D.length2 ppMatrix.Matrix)-1 do
                pwMatrix.[item, position] <- ppMatrix.[item, position]/pcv.[item]
        pwMatrix

    /// Calculate the score of the given sequence based on the Probabilities of the PositionWeightMatrix.
    let calculateSegmentScoreBy (pwMatrix:PositionWeightMatrix) (bioItems:BioArray.BioArray<#IBioItem>) =
        Array.fold (fun (position:int, value:float) (bios:#IBioItem) -> 
            position + 1, value * (pwMatrix.[bios, position])) (0, 1.) bioItems
        |> snd


open CompositeVector

module SiteSampler =

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSsWithBPV (motiveLength:int) (alphabet:#IBioItem[]) (source:BioArray.BioArray<#IBioItem>) (pcv:ProbabilityCompositeVector) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motiveLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motiveLength
                    let pwMatrix = PositionMatrix.createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                    segment
                    |> PositionMatrix.calculateSegmentScoreBy pwMatrix
                if tmp > highValue then loop (n + 1) tmp n
                else loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Checks whether downstream of given positions a higher InformationContent is present or not. 
    /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
    let getRightShiftedBestPWMSsWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    let tmp = Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motiveLength - 1 then position + 1 else position) tmp
                let unChosenArrays =
                    Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (PositionMatrix.getSegment motiveLength subSequence position) 
                         |> PositionMatrix.createPFMOf) unChosenArrays unChosenStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
    /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
    let getLeftShiftedBestPWMSsWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                let unChosenArrays =
                    Array.append sources.[0..randomSourceNumber.[randomSourceNumber.[n]]-1] sources.[randomSourceNumber.[randomSourceNumber.[n]]+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (PositionMatrix.getSegment motiveLength subSequence position) 
                         |> PositionMatrix.createPFMOf) unChosenArrays unChosenStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks the given Sequence for the existence of a conserved motive, by scoring each segment based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence.
    let findBestmotiveWithStartPosition (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =        
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) acc bestmotive =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..randomSourceNumber.[n]-1] acc.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> position)
                let unChosenArrays =
                    Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (PositionMatrix.getSegment motiveLength subSequence position) 
                         |> PositionMatrix.createPFMOf) unChosenArrays unChosenStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getPWMOfRandomStartsWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =    
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> List.toArray
            else
                let unChosenArrays =
                    Array.append (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        PositionMatrix.getRandomNumberInSequence motiveLength unChosen)
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (PositionMatrix.getSegment motiveLength subSequence position) 
                         |> PositionMatrix.createPFMOf) unChosenArrays randomStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                loop (n + 1) (getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix::acc)
        loop 0 []

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getmotivesWithBestInformationContentWithBPV (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStartsWithBPV motiveLength pseudoCount alphabet sources pcv
                            |> findBestmotiveWithStartPosition motiveLength pseudoCount alphabet sources pcv
                            |> getLeftShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv
                            |> getRightShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSs (motiveLength:int) (alphabet:#IBioItem[]) (pseudoCount:float) (source:BioArray.BioArray<#IBioItem>) (fcVector:FrequencyCompositeVector) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motiveLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motiveLength
                    let pcv =
                        increaseInPlaceFCVOf source fcVector
                        |> substractSegmentCountsFrom segment 
                        |> createNormalizedPCVOfFCV alphabet pseudoCount
                    let pwMatrix = PositionMatrix.createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                    segment
                    |> PositionMatrix.calculateSegmentScoreBy pwMatrix
                if tmp > highValue then loop (n + 1) tmp n
                else loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Checks whether downstream of given positions a higher InformationContent is present or not. 
    /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
    let getRightShiftedBestPWMSs (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    let tmp = Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motiveLength - 1 then position + 1 else position) tmp
                let unChosenArrays =
                    Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        PositionMatrix.getSegment motiveLength subSequence position
                        |> PositionMatrix.createPFMOf) unChosenArrays unChosenStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
    /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
    let getLeftShiftedBestPWMSs (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
            if n = sources.Length then
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                let unChosenArrays =
                    Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        PositionMatrix.getSegment motiveLength subSequence position
                        |> PositionMatrix.createPFMOf) unChosenArrays unChosenStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks the given sequence for the existence of a conserved motive, by scoring each segment based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated at each iteration if segments with higher scores are found until convergence.
    let getBestPWMSsWithStartPositions (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =        
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) acc bestmotive =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..randomSourceNumber.[n]-1] acc.[randomSourceNumber.[n]+1..]
                    |> Array.map (fun (_, position) -> position)
                let unChosenArrays =
                    Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        PositionMatrix.getSegment motiveLength subSequence position
                        |> PositionMatrix.createPFMOf) unChosenArrays unChosenStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                        acc.[randomSourceNumber.[n]] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getPWMOfRandomStarts (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) (acc:((float*int)[])) =
            if n = sources.Length then acc
            else
                let unChosenArrays =
                    Array.append (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        PositionMatrix.getRandomNumberInSequence motiveLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        PositionMatrix.getSegment motiveLength subSequence position
                        |> PositionMatrix.createPFMOf) unChosenArrays randomStartPositions
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                              acc)
        loop 0 (Array.zeroCreate sources.Length)

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getmotivesWithBestInformationContent (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStarts motiveLength pseudoCount alphabet sources
                            |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
                            |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                            |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getmotivesWithBestPWMSOfPPM (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix) =
        let randomSourceNumber = Array.shuffleFisherYates([|0..sources.Length-1|])
        let rec loop (n:int) acc:((float*int)[]) =
            if n = sources.Length then acc
            else
                let unChosenArrays =
                    Array.append (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        PositionMatrix.getRandomNumberInSequence motiveLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        createFCVWithout motiveLength position unchosenArray) 
                            unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                              acc)
        loop 0 (Array.zeroCreate sources.Length)

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getBestInformationContentOfPPM (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getmotivesWithBestPWMSOfPPM motiveLength pseudoCount alphabet sources positionProbabilityMatrix
                            |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
                            |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                            |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    let doSiteSamplingWithBPV (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        getPWMOfRandomStartsWithBPV motiveLength pseudoCount alphabet sources pcv
        |> findBestmotiveWithStartPosition motiveLength pseudoCount alphabet sources pcv
        |> getLeftShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv
        |> getRightShiftedBestPWMSsWithBPV motiveLength pseudoCount alphabet sources pcv

    let doSiteSampling (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        getPWMOfRandomStarts motiveLength pseudoCount alphabet sources
        |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources

    let doSiteSamplingWithPPM (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (ppM:PositionMatrix.PositionProbabilityMatrix) =
        getmotivesWithBestPWMSOfPPM motiveLength pseudoCount alphabet sources ppM
        |> getBestPWMSsWithStartPositions motiveLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motiveLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motiveLength pseudoCount alphabet sources

module motiveSampler =

    /// Contains the information of the probability and the start position(s) of the found motives(s). 
    type motiveIndex =
        {
            PWMS        : float
            Positions   : int list
        }

    /// Creates a motiveIndex based on the probability of the given segement(s) start position(s).
    let createmotiveIndex pwms pos =
        {
            motiveIndex.PWMS         = pwms
            motiveIndex.Positions    = pos
        }

    /// Calculates all possible segment combination of m segments, which do not overlap in the given width 
    /// and gives you the probability for it.
    let calculatePWMsForSegmentCombinations (cutOff:float) (width:int) (m:int) (set:list<float*int>) =
        let rec loop prob positions size set = 
            seq
                {
                    match size, set with
                    | n, x::xs ->
                        if n > 0 then
                            if PositionMatrix.ceckForDistance width (snd x::positions) then
                                if log2(fst x*prob) > cutOff then
                                    yield! loop (fst x*prob) (snd x::positions) (n - 1) xs
                        if n >= 0 then yield! loop prob positions n xs
                    | 0, [] -> yield createmotiveIndex (log2(prob)) positions
                    | _, [] -> () 
                }
        loop 1. [] m set
        |> List.ofSeq

    /// Normalizes the probabilities of all motiveMemories to the sum of all probabilities and picks one by random, 
    /// but those with a higher value have a higher chance to get picked.
    let rouletteWheelSelection (pick:float) (items:motiveIndex list) =
        let normalizedItems =
            let sum = List.sum (items |> List.map (fun item -> item.PWMS))
            items
            |> List.map (fun item -> item.PWMS/sum, item)
        let rec loop acc n =
            if acc <= pick && pick <= acc + fst normalizedItems.[n] then items.[n]
            else loop (acc + fst normalizedItems.[n]) (n + 1)
        loop 0. 0       

    /// Calculates the normalized segment score based on the given PositionWeightMatrix and 
    /// BackGroundProbabilityVecor of all segment combinations of a given sequence. 
    /// The amount of combinations is given by the amount of motives.
    let calculateNormalizedSegmentScores (cutOff:float) (motiveAmount:int) (motiveLength:int) (source:BioArray.BioArray<#IBioItem>) (pcv:ProbabilityCompositeVector) (pwMatrix:PositionMatrix.PositionWeightMatrix) =
        let segments =
            let rec loop n acc =
                if n + motiveLength = source.Length+1 then List.rev acc
                else
                    let tmp =
                        source
                        |> Array.skip n
                        |> (fun items -> Array.take motiveLength items, n)
                    loop (n+1) (tmp::acc)
            loop 0 []
        let segmentScores =
            segments
            |> List.map (fun segment -> 
                PositionMatrix.calculateSegmentScoreBy pwMatrix (fst segment), (snd segment))
        let backGroundScores =
            segments
            |> List.map (fun segment -> calculateSegmentScoreBy pcv (fst segment))
            |> List.map (fun segmentScore -> createmotiveIndex segmentScore [])
        let tmp =
            let rec loop n acc =
                if n > motiveAmount then List.rev acc
                else loop (n+1) (calculatePWMsForSegmentCombinations cutOff motiveLength n segmentScores::acc)
            loop 1 [backGroundScores]
        tmp
        |> List.concat

    /// Checks the given sequence for the existence of conserved motive(s), by scoring each segment and segment combination based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence. 
    let findBestmotivePositionsWithStartPositionByPCV (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (motiveMem:motiveIndex[]) =
        let rec loop (n:int) acc bestmotive =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> item.Positions)) = (bestmotive |> Array.map (fun item -> item.Positions)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = PositionMatrix.createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motiveAmount motiveLength sources.[n] pcv positionWeightMatrix
                    |> List.sortByDescending (fun item -> item.PWMS)
                    |> List.head
                loop 
                    (n + 1) 
                    (if tmp.PWMS > acc.[n].PWMS then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy motiveMem) (Array.copy motiveMem)


    /// Checks the given Sequence for the existence of conserved motive(s), by scoring each segment and segment combination based on the given start positions.
    /// The nex segments for the nex PositionWeightMatrix are picked by chance after calculating the segment scores for each possible combination 
    /// but those with higher scores have a higher chance to be picked.
    let findBestmotivePositionsWithStartPositionsByPCV (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (motiveMem:motiveIndex[]) =
        let rnd = new System.Random()
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> Array.ofList
            else
                let unChosenStartPositions =
                    Array.append motiveMem.[0..n-1] motiveMem.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = PositionMatrix.createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motiveAmount motiveLength sources.[n] pcv positionWeightMatrix
                    |> rouletteWheelSelection (rnd.NextDouble())
                loop (n+1) (tmp::acc)
        loop 0 []

    ///Repeats the motive sampler until convergence or until a given number of repetitions is done.
    let findBestInormationContentContainingmotivesWithPCV (numberOfRepetitions:int) (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        let rec loop (n:int) (acc:motiveIndex[]) (bestPWMS:motiveIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            SiteSampler.getPWMOfRandomStartsWithBPV motiveLength pseudoCount alphabet sources pcv
                            |> Array.map (fun (prob, position) -> createmotiveIndex prob [position])
                            |> findBestmotivePositionsWithStartPositionsByPCV motiveAmount motiveLength pseudoCount cutOff alphabet sources pcv
                            |> findBestmotivePositionsWithStartPositionByPCV motiveAmount motiveLength pseudoCount cutOff alphabet sources pcv
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createmotiveIndex 0. []|]

    /// Checks the given Sequence for the existence of conserved motive(s), by scoring each segment and segment combination based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence. 
    let findBestmotiveIndicesWithStartPositions (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motiveMem:motiveIndex[]) =
        let rec loop (n:int) acc bestmotive =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> item.Positions)) = (bestmotive |> Array.map (fun item -> item.Positions)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions
                        |> Array.map (fun position -> 
                            createFCVWithout motiveLength position array)
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> increaseInPlaceFCVOf sources.[n]
                    |> createNormalizedPCVOfFCV alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = PositionMatrix.createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motiveAmount motiveLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> List.sortByDescending (fun item -> item.PWMS)
                    |> List.head
                loop 
                    (n + 1) 
                    (if tmp.PWMS > acc.[n].PWMS then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestmotive
        loop 0 (Array.copy motiveMem) (Array.copy motiveMem)


    /// Checks the given Sequence for the existence of conserved motive(s), by scoring each segment and segment combination based on the given start positions.
    /// The nex segments for the nex PositionWeightMatrix are picked by chance after calculating the segment scores for each possible combination 
    /// but those with higher scores have a higher chance to be picked.
    let findBestmotiveIndicesByWithStartPositions (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motiveMem:motiveIndex[]) =
        let rnd = new System.Random()
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> Array.ofList
            else
                let unChosenStartPositions =
                    Array.append motiveMem.[0..n-1] motiveMem.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions
                        |> Array.map (fun position -> 
                            createFCVWithout motiveLength position array)
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> increaseInPlaceFCVOf sources.[n]
                    |> createNormalizedPCVOfFCV alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOf)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                    |> PositionMatrix.createPPMOf
                    |> PositionMatrix.normalizePPM (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = PositionMatrix.createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motiveAmount motiveLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> rouletteWheelSelection (rnd.NextDouble())
                loop (n+1) (tmp::acc)
        loop 0 []

    ///Repeats the motive sampler until convergence or until a given number of repetitions is done.
    let getmotivesWithBestInformationContents (numberOfRepetitions:int) (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) (acc:motiveIndex[]) (bestPWMS:motiveIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            SiteSampler.getPWMOfRandomStarts motiveLength pseudoCount alphabet sources
                            |> Array.map (fun (prob, position) -> createmotiveIndex prob [position])
                            |> findBestmotiveIndicesByWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources
                            |> findBestmotiveIndicesWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createmotiveIndex 0. []|]

    ///Repeats the motive sampler until convergence or until a given number of repetitions is done.
    let getBestPWMSsOfPPM (numberOfRepetitions:int) (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix) =
        let rec loop (n:int) (acc:motiveIndex[]) (bestPWMS:motiveIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            SiteSampler.getmotivesWithBestPWMSOfPPM motiveLength pseudoCount alphabet sources positionProbabilityMatrix
                            |> Array.map (fun (prob, position) -> createmotiveIndex prob [position])
                            |> findBestmotiveIndicesByWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources
                            |> findBestmotiveIndicesWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createmotiveIndex 0. []|]

    let domotiveSamplingWithPPM (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix) =
        SiteSampler.getmotivesWithBestPWMSOfPPM motiveLength pseudoCount alphabet sources positionProbabilityMatrix
        |> Array.map (fun (prob, position) -> createmotiveIndex prob [position])
        |> findBestmotiveIndicesByWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources
        |> findBestmotiveIndicesWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources

    let domotiveSampling (motiveAmount:int) (motiveLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        SiteSampler.getPWMOfRandomStarts motiveLength pseudoCount alphabet sources
        |> Array.map (fun (prob, position) -> createmotiveIndex prob [position])
        |> findBestmotiveIndicesByWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources
        |> findBestmotiveIndicesWithStartPositions motiveAmount motiveLength pseudoCount cutOff alphabet sources

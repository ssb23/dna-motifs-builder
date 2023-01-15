from flask import Flask, render_template, request, redirect, url_for, make_response
import generateMotifs
#import pdfkit
import asyncio

app = Flask(__name__)



@app.route('/', methods=['GET', 'POST'])
async def validate():
    form = {'gcContentMinPercentage': 20, 
            'gcContentMaxPercentage': 60, 
            'payloadSize': 5, 
            'payloadNum': 10,
            'keySize': 2,
            'keyNum': 2,
            'maxHomopolymer':1,
            'maxHairpin': 2,
            'loopMin': 2,
            'loopMax': 2,
            'keyGc': False, 
            'homVisible': False,
            'hairpinVisible': False,
            'gcVisible': False,
            }

    if request.method == 'POST':

        if "analyseSeqSubmission" in request.form:
            # Form being submitted; grab data from form.
            gcContentMinPercentage = 20
            gcContentMaxPercentage = 60
            payloadSize = int(request.form['payloadSize'])
            payloadNum = int(request.form['payloadNum'])
            keyNum = int(request.form['keyNum'])
            keySize = int(request.form['keySize'])
            maxHomopolymer = 1
            maxHairpin = 2
            loopMin = 2
            loopMax = 2
            constraints = set()
            if request.form['homVisible'] == 'True':
                constraints.add('hom')
                maxHomopolymer = int(request.form['maxHomopolymer'])
            if request.form['hairpinVisible'] == 'True':
                constraints.add('hairpin')
                maxHairpin = int(request.form['maxHairpin'])
                loopMin = int(request.form['loopMin'])
                loopMax = int(request.form['loopMax'])
            if request.form['gcVisible'] == 'True':
                constraints.add('motifGcContent')
                constraints.add('keyGcContent')
                gcContentMinPercentage = int(request.form['gcContentMinPercentage'])
                gcContentMaxPercentage = int(request.form['gcContentMaxPercentage'])

            form = {'gcContentMinPercentage': gcContentMinPercentage, 
                    'gcContentMaxPercentage': gcContentMaxPercentage, 
                    'payloadSize': payloadSize, 
                    'payloadNum': payloadNum,
                    'keySize': keySize,
                    'keyNum': keyNum,
                    'maxHomopolymer': maxHomopolymer,
                    'maxHairpin': maxHairpin,
                    'loopMin': loopMin,
                    'loopMax': loopMax,
                    'keyGc': False, 
                    'homVisible': False,
                    'hairpinVisible': False,
                    'gcVisible': False, 
                    # 'keyGc': request.form['keyGc'], 
                    # 'homVisible': request.form['homVisible'],
                    # 'hairpinVisible': request.form['hairpinVisible'],
                    # 'gcVisible': request.form['gcVisible'],
                    }

            # Default values
            motifs = ""
            keys = ""
            isValid = False

            blob = await asyncio.gather(generateMotifs.generateMotifs(constraints, payloadNum, payloadSize, keyNum, keySize, maxHairpin, gcContentMinPercentage, gcContentMaxPercentage, maxHomopolymer, loopMin, loopMax))
            stringBlob = str(blob)
            if stringBlob[len(stringBlob)-4] == 's':
                setThing = stringBlob[2:len(stringBlob)-9]
                isValid  = False
            else:
                setThing = stringBlob[2:len(stringBlob)-8]
                isValid  = True

            keys, payloads, motifs = setThing.split('}, {')
            keys = keys[2:len(keys) - 1]
            payloads = payloads[1:len(payloads) - 1]
            motifs = motifs[1:len(motifs) - 2]

            finalIsValid = isValid

            finalKeys = ""
            for k in keys.split("', '"):
                finalKeys += k + "    "

            finalPayloads = ""
            for p in payloads.split("', '"):
                finalPayloads += p + "    "

            finalMotifs = ""
            for m in motifs.split("', '"):
                finalMotifs += m + "    "
            # Render the page
            return render_template('index.html', payloads=finalPayloads, motifs=finalMotifs, keys=finalKeys, form=form, isValid=isValid)
        
        return render_template('index.html', payloads="", motifs="", keys="", form=form, isValid=False)

    return render_template('index.html', payloads="", motifs="", keys="", form=form, isValid=False)


if __name__ == "__main__":
   app.run(debug=True)



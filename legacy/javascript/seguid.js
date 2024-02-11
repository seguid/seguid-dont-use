function makeComplementTable() {
    const dict = {
        "A": "T", "C": "G", "G": "C",
        "T": "A", "M": "K", "R": "Y",
        "W": "W", "S": "S", "Y": "R",
        "K": "M", "V": "B", "H": "D",
        "D": "H", "B": "V", "X": "X",
        "N": "N", "U": "A",
    };
    const keys = Object.keys(dict).join("") + Object.keys(dict).join("").toLowerCase();
    const values = Object.values(dict).join("") + Object.values(dict).join("").toLowerCase();
    const table = {};

    for (let i = 0; i < keys.length; i++) {
        table[keys.charCodeAt(i)] = values[i];
    }

    return table;
}

const complement_table = makeComplementTable();

const purify = (s) => {
    return s.replace(/[^a-zA-Z]/g, '').toUpperCase();
}

const complement = (s) => {
    const chars = s.split('');
    for (let i = 0; i < chars.length; i++) {
        const charCode = chars[i].charCodeAt(0);
        if (complement_table[charCode]) {
            chars[i] = complement_table[charCode];
        }
    }
    return chars.join('');
}

const reverse = (s) => {
    return s.split('').reverse().join('');
}

const reverse_complement = (s) => {
    return complement(reverse(s));
}

const SEGUID = async (seq) => {
    const encoder = new TextEncoder();
    const data = encoder.encode(seq.toUpperCase());
    try {
        const hashBuffer = await crypto.subtle.digest('SHA-1', data);
        const hashArray = Array.from(new Uint8Array(hashBuffer));
        const hashString = hashArray.map(b => String.fromCharCode(b)).join('');
        const encodedHash = btoa(hashString);
        return encodedHash.replace(/\n|=/g, '').replace(/\+/g, '-').replace(/\//g, '_');
    } catch (error) {
        console.error('Error calculating SHA-1 hash:', error);
        return null;
    }
};

const minRotation = (s) => {
    let N = s.length;
    s += s;
    let a = 0, b = 0;
    while (b < N) {
        for (let i = 0; i < N - a; i++) {
            let sai = s[a + i];
            let sbi = s[b + i];
            if (sai < sbi || a + i === b) {
                if (i) {
                    b += i - 1;
                }
                break;
            }
            if (sai > sbi) {
                a = b;
                break;
            }
        }
        b += 1;
    }
    return s.slice(a, a + N);
}

const lSEGUID = (s) => {
    s = s.toUpperCase();
    const rc = reverse_complement(s);
    if (rc < s)
        s = rc;
    return SEGUID(s);
}
const cSEGUID = (s) => {
    s = s.toUpperCase();
    const rc = reverse_complement(s);
    const rot_s = minRotation(s);
    const rot_rc = minRotation(rc);
    return SEGUID(rot_s < rot_rc ? rot_s : rot_rc);
}
